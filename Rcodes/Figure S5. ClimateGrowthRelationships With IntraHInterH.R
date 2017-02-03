rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
load(file.path(workPath, "data", selectionMethod, "fullClimateModels.RData"))

competitionModel <- c("indiH")
for(i in competitionModel){
  allcoeffTable <- get(paste(i, "FixedCoeff", sep = ""))
  for(j in 1:length(allcoeffTable)){
    modelName <- unlist(strsplit(names(allcoeffTable)[j], "_", fixed = T))
    coeffIndi <- cbind(data.table(competitionModel = i,
                                  Model = names(allcoeffTable)[j],
                                  Species = modelName[1],
                                  Climate = modelName[2]), 
                       allcoeffTable[[j]])
    if(i == "indiH" & j == 1){
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
for(i in c("IntraH", "InterH")){
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

interactionTable <- OverallResults[rn %in% c("Climatectd:logIntraHctd", "Climatectd:logInterHctd"),][
  ,.(competitionModel, rn, Species, Climate, interactEff = Value, 
     interactEff_SE = Std.Error, Pvalue = `p-value`)]
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
alloutput1[,lineTransp:=1]
alloutput1[Pvalue > 0.05, lineTransp:=2]

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
climateWithCompTable[Species == "Minor species", Species:="Minor species group"]
climateWithCompTable[Species == "All species", Species:="All trees"]
majorspecies <- c("All trees", majorspecies[2:4], "Minor species group")


climateWithCompTable[,':='(Species = factor(Species, levels = majorspecies))]

climateWithCompTable[Climate %in% c("ATA", "GSTA", "NONGSTA"), ClimateName:="Temperature"]
climateWithCompTable[Climate %in% c("ACMIA", "GSCMIA", "NONGSCMIA"), ClimateName:="CMI"]
climateWithCompTable[Climate %in% c("ACO2A", "GSCO2A", "NONGSCO2A"), ClimateName:="CO2"]
climateWithCompTable[, ClimateName:=factor(ClimateName, levels = c("Temperature", "CMI", "CO2"))]



ZeroLines <- rbind(data.table(SeasonComp = c("Annual anomaly", 
                                             "Growing season anomaly", 
                                             "Non-growing season anomaly"),
                              ClimateName="Temperature"),
                   data.table(SeasonComp = c("Annual anomaly", 
                                             "Growing season anomaly"),
                              ClimateName="CMI"),
                   data.table(SeasonComp = c("Annual anomaly"),
                              ClimateName="CO2"))

ZeroLines[,':='(SeasonComp = factor(SeasonComp, levels = c("Annual anomaly", 
                                                           "Growing season anomaly", 
                                                           "Non-growing season anomaly")),
                ClimateName = factor(ClimateName, levels = c("Temperature", "CMI", "CO2")))]
ZeroLines[,':='(x = -Inf, xend = Inf, y = 0, yend = 0)]

HbreakLines <- rbind(ZeroLines[,.(SeasonComp, ClimateName)][,':='(x = 6, xend = 6, y = -Inf, yend = Inf)],
                     ZeroLines[,.(SeasonComp, ClimateName)][,':='(x = -Inf, xend = -Inf, y = -Inf, yend = Inf)])

Yaxislines <- data.table(SeasonComp = factor(c("Annual anomaly"),
                                             levels = c("Annual anomaly", 
                                                        "Growing season anomaly", 
                                                        "Non-growing season anomaly")),
                         x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

Xaxislines <- data.table(x = Inf, xend = -Inf, y = -Inf, yend = -Inf)

labeltexts <- data.table(ClimateName = factor(c("Temperature", "CMI", "CO2"),
                                              levels = c("Temperature", "CMI", 
                                                         "CO2")), 
                         SeasonComp = factor(c("Annual anomaly"),
                                             levels = c("Annual anomaly", 
                                                        "Growing season anomaly", 
                                                        "Non-growing season anomaly")),
                         labels = letters[1:3],
                         x = -Inf, y = Inf)

newlabels1 <- list("Temperature" = "Temperature anomaly effect",
                   "CMI" = "CMI anomaly effect",
                   'CO2'= expression(paste("C", O[2], " anomaly effect")))
newlabels2 <- list("Annual anomaly" = "Annual",
                   "Growing season anomaly" = "Growing season",
                   "Non-growing season anomaly" = "Non-growing season")


figure_labeller <- function(variable,value){
  if(variable == "ClimateName"){
    return(newlabels1[value])
  } else if (variable == "SeasonComp"){
    return(newlabels2[value])
  } 
}
maxTemperatureEff <- 1.2
IntraHInterHTexts <- data.table(ClimateName = factor("Temperature", levels = c("Temperature", "CMI", 
                                                                               "CO2")),
                                x = c(3, 3+6),
                                y = maxTemperatureEff,
                                texts = c("IntraH", "InterH"))

climateWithCompTable[EffectMin>=EffectMax,
                     ':='(EffectStart = EffectMin+1.98*EffectMin_SE,
                          EffectEnd = EffectMax-1.98*EffectMax_SE)]
climateWithCompTable[EffectMin<EffectMax,
                     ':='(EffectStart = EffectMin-1.98*EffectMin_SE,
                          EffectEnd = EffectMax+1.98*EffectMax_SE)]

FigureB <- ggplot(data = climateWithCompTable[lineTransp==1,], aes(x = xaxis, y = mainEffect))+
  facet_grid(ClimateName~SeasonComp,
             scales = "free_y", switch = "y",
             labeller = figure_labeller, drop = FALSE)+
  geom_segment(data = ZeroLines, aes(x = x, xend = xend, 
                                     y = y, yend = yend), linetype = 2, col = "gray", size = 1)+
  geom_segment(data = HbreakLines[x != -Inf,], aes(x = x, xend = xend, 
                                     y = y, yend = yend), linetype = 1, col = "gray", size = 1)+
  # geom_errorbar(aes(ymin = EffectMin-1.98*EffectMin_SE, ymax = EffectMin+1.98*EffectMin_SE), 
  #               width = 0.4, col = "gray", size = 1)+
  # 
  geom_segment(aes(col = Species, y = EffectStart,
                   yend = EffectEnd, 
                   x = xaxis, xend = xaxis), 
               arrow = arrow(length = unit(0.05, "npc")), size = 1)+
  geom_errorbar(data = climateWithCompTable[lineTransp==2,],
                aes(ymin = pmin(EffectStart, EffectEnd),
                    ymax = pmax(EffectStart, EffectEnd),
                    col = Species,), 
                               width = 0.4, size = 1)+

  # geom_errorbar(data = climateWithCompTable[interactEff == 0 & mainEffect != 0, ],
  #               aes(ymin = mainEffect-1.98*mainEffect_SE, ymax = mainEffect+1.98*mainEffect_SE),
  #               col = "gray", size = 1, width = 0.2)+
  scale_x_continuous(name = "a", limits = c(0, 12))+
  geom_point(data = climateWithCompTable[interactEff == 0 & mainEffect != 0, ], 
             aes(y = mainEffect, col = Species),
             size = 2, show.legend = FALSE)+
  geom_segment(data = Yaxislines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1.5, col = "black")+
  geom_segment(data = Xaxislines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1, col = "black")+
  geom_text(data = labeltexts, aes(x = x, y = y, label = labels),
            vjust = 1.5, hjust = -1.2, size = 7)+
  geom_text(data = IntraHInterHTexts, aes(x = x, y = y, label = texts), size = 3)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85, 0.50),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.background = element_rect(colour = "black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 15, vjust = -2))
workPath <- "~/GitHub/Climate_Growth"

ggsave(file.path(workPath, "TablesFigures",
                 "Figure S5. Climate associations with IntraHInterH.png"), FigureB,
       width = 10, height = 7)


