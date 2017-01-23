rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "ClimateModels_AllCensus_PositiveGrowth.RData"))

competitionModel <- c("allH", "indiH")
for(i in competitionModel){
  allcoeffTable <- get(paste(i, "FixedCoeff", sep = ""))
  for(j in 1:length(allcoeffTable)){
    modelName <- unlist(strsplit(names(allcoeffTable)[j], "_", fixed = T))
    coeffIndi <- cbind(data.table(competitionModel = i,
                                  Model = names(allcoeffTable)[j],
                                  Species = modelName[1],
                                  Climate = modelName[2]), 
                       allcoeffTable[[j]])
    if(i == "allH" & j == 1){
      allCoeff <- coeffIndi
    } else {
      allCoeff <- rbind(allCoeff, coeffIndi)
    }
  }
}
for(i in 1:nrow(allCoeff)){
  allCoeff$rn[i] <- paste(sort(unlist(strsplit(allCoeff$rn[i], ":", fixed = TRUE))), collapse = ":")
}


OverallResults <- data.table::copy(allCoeff)
temperature <- c("ATA", "GSTA", "NONGSTA")
CMI <- c("ACMIA", "GSCMIA")
CO2 <- c("ACO2A")
longcol <- c(temperature, CMI, CO2)
majorspecies <- studySpecies
for(indispecies in majorspecies){
  speciesData <- analysesData[Species == indispecies,]
  climateWithCompetitionTable <- rbind(data.table(expand.grid(Species = indispecies,
                                                              competitionModel = competitionModel,
                                                          Climate = as.character(longcol), 
                                                          CompetitionName = "H",
                                                          CompetitionMin = min(speciesData$H),
                                                          CompetitionMax = max(speciesData$H),
                                                          stringsAsFactors = FALSE)),
                                   data.table(expand.grid(Species = indispecies,
                                                          competitionModel = competitionModel,
                                                          Climate = as.character(longcol), 
                                                          CompetitionName = "IntraH",
                                                          CompetitionMin = min(speciesData$IntraH),
                                                          CompetitionMax = max(speciesData$IntraH),
                                                          stringsAsFactors = FALSE)),
                                   data.table(expand.grid(Species = indispecies,
                                                          competitionModel = competitionModel,
                                                          Climate = as.character(longcol), 
                                                          CompetitionName = "InterH",
                                                          CompetitionMin = min(speciesData$InterH),
                                                          CompetitionMax = max(speciesData$InterH),
                                                          stringsAsFactors = FALSE))) %>% data.table
  climateWithCompetitionTable[CompetitionName == "H",
                              ':='(CompetitionMinctd = log(CompetitionMin+1)-mean(log(speciesData$H+1)),
                                   CompetitionMaxctd = log(CompetitionMax+1)-mean(log(speciesData$H+1)))]
  climateWithCompetitionTable[CompetitionName == "IntraH",
                          ':='(CompetitionMinctd = log(CompetitionMin+1)-mean(log(speciesData$IntraH+1)),
                               CompetitionMaxctd = log(CompetitionMax+1)-mean(log(speciesData$IntraH+1)))]
  climateWithCompetitionTable[CompetitionName == "InterH",
                          ':='(CompetitionMinctd = log(CompetitionMin+1)-mean(log(speciesData$InterH+1)),
                               CompetitionMaxctd = log(CompetitionMax+1)-mean(log(speciesData$InterH+1)))]
  if(indispecies == "All species"){
    alloutput <- climateWithCompetitionTable
  } else {
    alloutput <- rbind(alloutput, climateWithCompetitionTable)
  }
}
m <- 1
for(i in c("H", "IntraH", "InterH")){
  alloutput[CompetitionName == i & Species == majorspecies[1], xaxis := m]
  alloutput[CompetitionName == i & Species == majorspecies[2], xaxis := m+1]
  alloutput[CompetitionName == i & Species == majorspecies[3], xaxis := m+2]
  alloutput[CompetitionName == i & Species == majorspecies[4], xaxis := m+3]
  alloutput[CompetitionName == i & Species == majorspecies[5], xaxis := m+4]
  m <- m+6
}


maineffectTable <- OverallResults[rn == "Climatectd",][, .(competitionModel,
                                                           Species, Climate, mainEffect = Value,
                                                           mainEffect_SE = Std.Error)]

interactionTable <- OverallResults[rn %in% c("Climatectd:logHctd", "Climatectd:logIntraHctd", "Climatectd:logInterHctd"),][
  ,.(competitionModel, rn, Species, Climate, interactEff = Value, 
     interactEff_SE = Std.Error)]
interactionTable[rn == "Climatectd:logHctd", CompetitionName:="H"]
interactionTable[rn == "Climatectd:logIntraHctd", CompetitionName:="IntraH"]
interactionTable[rn == "Climatectd:logInterHctd", CompetitionName:="InterH"]
interactionTable[,rn:=NULL]


alloutput1 <- dplyr::left_join(alloutput, maineffectTable, by = c("competitionModel", "Species", "Climate")) %>%
  data.table
alloutput1[is.na(mainEffect), ':='(mainEffect = 0, mainEffect_SE = 0)]

alloutput1 <- dplyr::left_join(alloutput1, interactionTable, by = c("competitionModel", "Species", "Climate", "CompetitionName")) %>%
  data.table
alloutput1 <- alloutput1[is.na(interactEff),':='(interactEff = 0, 
                                                 interactEff_SE = 0)]
alloutput1[,':='(EffectMin = mainEffect+CompetitionMinctd*interactEff, 
                 EffectMin_SE = sqrt(mainEffect_SE^2+interactEff_SE^2),
                 EffectMax = mainEffect+CompetitionMaxctd*interactEff, 
                 EffectMax_SE = sqrt(mainEffect_SE^2+interactEff_SE^2))]

climateWithCompTable <- data.table::copy(alloutput1)
climateWithCompTable[Climate %in% c("ATA", "ACMIA", "ACO2A"),
                     SeasonComp:="Annual anomaly"]
climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A"),
                     SeasonComp:="Growing season anomaly"]
climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A"),
                     SeasonComp:="Non-growing season anomaly"]

climateWithCompTable[, SeasonComp:=factor(SeasonComp, 
                                          levels = c("Annual anomaly", 
                                                     "Growing season anomaly", 
                                                     "Non-growing season anomaly"))]
climateWithCompTable[Species == "Other species", Species:="Minor species"]
majorspecies <- c(majorspecies[1:4], "Minor species")

climateWithCompTable[,':='(Species = factor(Species, levels = majorspecies))]

climateWithCompTable[Climate %in% c("ATA", "GSTA", "NONGSTA"), ClimateName:="Temperature"]
climateWithCompTable[Climate %in% c("ACMIA", "GSCMIA", "NONGSCMIA"), ClimateName:="CMI"]
climateWithCompTable[Climate %in% c("ACO2A", "GSCO2A", "NONGSCO2A"), ClimateName:="CO2"]
climateWithCompTable[, ClimateName:=factor(ClimateName, levels = c("Temperature", "CMI", "CO2"))]

climateWithCompTable <- climateWithCompTable[Climate %in% c("ATA", "ACMIA", "ACO2A"),]


# startPoints <- climateWithCompTable[xscale %in% c(1, 4, 7, 10),][, ':='(x = xscale+1)]
# endPoints <- climateWithCompTable[xscale %in% c(3, 6, 9, 12), ][, .(Species, Climate, CompetitionName, 
#                                                                     xend = xscale-1, 
#                                                                     Effectend = Effect, Effectend_SE = Effect_SE)]
# segmentPoints <- setkey(startPoints, Species, Climate, CompetitionName)[setkey(endPoints, Species, Climate, CompetitionName), 
#                                                                         nomatch = 0]
# segmentPoints[,':='(Species = factor(Species, levels = majorspecies))]
# mainEffect <- climateWithCompTable[xscale %in% c(2, 5, 8, 11) & mainEffect != 0,][,':='(x = xscale)]
# mainEffect[,':='(Species = factor(Species, levels = majorspecies))]

newlabels1 <- list("Temperature" = "Temperature effect",
                   "CMI" = "CMI effect",
                   'CO2'= expression(paste("C", O[2], " effect")))
newlabels2 <- list("WholeIntraH" = "IntraH",
                   "WholeInterH" = "InterH",
                   "GSIntraH" = "IntraH",
                   "GSInterH" = "InterH",
                   "NGSIntraH" = "IntraH",
                   "NGSInterH" = "InterH")


figure_labeller <- function(variable,value){
  if(variable == "ClimateName"){
    return(newlabels1[value])
  } else if (variable == "SeasonComp"){
    return(newlabels2[value])
  } 
}
majorYpanelbreaklines <- data.table(SeasonComp = factor(c("WholeIntraH"),
                                                        levels = c("WholeIntraH", "WholeInterH", 
                                                                   "GSIntraH", "GSInterH", 
                                                                   "NGSIntraH", "NGSInterH")),
                                    x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

majorXpanelbreaklines <- data.table(x = Inf, xend = -Inf, y = -Inf, yend = -Inf)
a <- unique(segmentPoints[,.(SeasonComp, ClimateName)],
            by = c("SeasonComp", "ClimateName"))
a[,':='(x = -1, xend = Inf, Effect = 0, yend = 0)]
labeltexts <- data.table(ClimateName = factor(c("Temperature", "CMI", "CO2"),
                                              levels = c("Temperature", "CMI", 
                                                         "CO2")), 
                         SeasonComp = factor(c("WholeIntraH"),
                                             levels = c("WholeIntraH", "WholeInterH", 
                                                        "GSIntraH", "GSInterH", 
                                                        "NGSIntraH", "NGSInterH")),
                         labels = letters[1:3],
                         x = 0.5, y = Inf)
seasontexts <- data.table(ClimateName = factor(c("Temperature"),
                                               levels = c("Temperature", "CMI", 
                                                          "CO2")), 
                          SeasonComp = factor(c("WholeIntraH", "GSIntraH", "NGSIntraH"),
                                              levels = c("WholeIntraH", "WholeInterH", 
                                                         "GSIntraH", "GSInterH", 
                                                         "NGSIntraH", "NGSInterH")),
                          labels = c("Annual anomaly", "Growing season anomaly",
                                     "Non-growing season anomaly"),
                          x = 2, y = Inf)

FigureB <- ggplot(data = climateWithCompTable[interactEff != 0, ], aes(x = xaxis, y = mainEffect))+
  facet_grid(ClimateName~CompetitionName,
             scales = "free_y", switch = "both")+
             # labeller = figure_labeller, drop = FALSE)+
  geom_point(data = climateWithCompTable[interactEff == 0 & mainEffect != 0, ], 
             aes(y = mainEffect, col = Species),
             size = 0.5)+
  geom_segment(aes(col = Species, y = EffectMin, yend = EffectMax, 
                   x = xaxis, xend = xaxis), 
               arrow = arrow(length = unit(0.05, "npc")), size = 1)+
  geom_segment(data = a, aes(x = x, xend = xend, 
                             y = Effect, yend = yend), linetype = 2, col = "gray", size = 1)+
  geom_errorbar(aes(group = Species, ymin = Effect-1.98*Effect_SE, ymax = Effect+1.98*Effect_SE), 
                width = 1.5, col = "gray", size = 1)+
  geom_errorbar(aes(group = Species, ymin = Effectend-1.98*Effectend_SE, ymax = Effectend+1.98*Effectend_SE), 
                width = 1.5, col = "gray", size = 1)+
  geom_errorbar(data = mainEffect[interactEff == 0, ], aes(x = x, ymin = mainEffect-mainEffect_SE, 
                                                           ymax = mainEffect+mainEffect_SE, group = Species),
                col = "gray", width = 1.5, size = 1)+
  geom_segment(aes(group = Species, col = Species, xend = xend, yend = Effectend), 
               arrow = arrow(length = unit(0.1, "npc")), size = 1)+
  geom_segment(data = majorYpanelbreaklines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1.5, col = "black")+
  geom_segment(data = majorXpanelbreaklines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1, col = "black")+
  geom_point(data = mainEffect[interactEff == 0, ], aes(x = x, y = mainEffect, group = Species, col = Species),
             show.legend = FALSE)+
  geom_text(data = labeltexts, aes(x = x, y = y, label = labels),
            vjust = 1, size = 7)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        # axis.title.x = element_text(colour = "black", size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85, 0.45),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.background = element_rect(colour = "black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 15, vjust = -2))


ggsave(file.path(workPath, "TablesFigures",
                 "Figure 4. Climate associations with tree growth.png"), FigureB,
       width = 10, height = 7)

library(png);library(grid); library(gridExtra)
img <- readPNG(file.path(workPath, "TablesFigures",
                         "Figure 4. Climate associations with tree growth.png"))
g <- rasterGrob(img)
g$width <- unit(1, "npc")
g$height <- unit(1, "npc")
start <- 0.24
dis <- 0.3
figureC <- ggplot(data = data.frame(x = c(0, 1), y = c(0, 1)),
                  aes(x = x, y = y))+
  annotation_custom(g, xmin = 0, xmax = 1, ymin = 0, ymax = 1)+
  geom_text(data = data.frame(x = c(start, start+dis, start+2*dis), y = 1, 
                              texts = c("Annual anomaly", "Growing season anomaly", "Non-growing season anomaly")),
            aes(x = x, y = y, label = texts), size = 5)+
  
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(-10,-25,-25,-25),
        panel.margin = margin(0,0,0,0))

ggsave(file.path(workPath, "TablesFigures",
                 "Figure 4. Climate associations with tree growth.png"), figureC,
       width = 10, height = 7)


