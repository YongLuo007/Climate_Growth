rm(list = ls())
library(data.table); library(ggplot2, lib.loc = "~/GitHub/Climate_Growth/RequiredRPackages");
# library(SpaDES, lib.loc = "~/GitHub/Climate_Growth/RequiredRPackages")
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data", selectionMethod, "fullClimateModels_BA.RData"))
useLogQunatile95 <- TRUE
competitionModel <- c("allH")
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
  if(useLogQunatile95){
    H95quantilelower <- exp(as.numeric(quantile(log(speciesData$MidH), probs = 0.025)))
    H95quantileupper <- exp(as.numeric(quantile(log(speciesData$MidH), probs = 0.975)))
    climateWithCompetitionTable <- rbind(data.table(expand.grid(Species = indispecies,
                                                                competitionModel = competitionModel,
                                                                Climate = as.character(longcol), 
                                                                CompetitionName = "H",
                                                                CompetitionMin = H95quantilelower,
                                                                CompetitionMax = H95quantileupper,
                                                                stringsAsFactors = FALSE))) %>% data.table
  } else {
    climateWithCompetitionTable <- rbind(data.table(expand.grid(Species = indispecies,
                                                                competitionModel = competitionModel,
                                                                Climate = as.character(longcol), 
                                                                CompetitionName = "H",
                                                                CompetitionMin = min(speciesData$MidH),
                                                                CompetitionMax = max(speciesData$MidH),
                                                                stringsAsFactors = FALSE))) %>% data.table
  }
  
  climateWithCompetitionTable[CompetitionName == "H",
                              ':='(CompetitionMinctd = log(CompetitionMin)-mean(log(speciesData$MidH)),
                                   CompetitionMaxctd = log(CompetitionMax)-mean(log(speciesData$MidH)))]
  if(indispecies == "All species"){
    alloutput <- climateWithCompetitionTable
  } else {
    alloutput <- rbind(alloutput, climateWithCompetitionTable)
  }
}
m <- 1
for(i in c("H")){
  alloutput[CompetitionName == i & Species == majorspecies[1], xaxis := m]
  alloutput[CompetitionName == i & Species == majorspecies[2], xaxis := m+1]
  alloutput[CompetitionName == i & Species == majorspecies[3], xaxis := m+2]
  alloutput[CompetitionName == i & Species == majorspecies[4], xaxis := m+3]
  alloutput[CompetitionName == i & Species == majorspecies[5], xaxis := m+4]
  m <- m+6
}


maineffectTable <- OverallResults[rn == "Climatectd",][, .(competitionModel,
                                                           Species, Climate, mainEffect = Value,
                                                           mainEffect_SE = Std.Error,
                                                           MainPvalue = `p-value`)]

interactionTable <- OverallResults[rn %in% c("Climatectd:logHctd"),][
  ,.(competitionModel, rn, Species, Climate, interactEff = Value, 
     interactEff_SE = Std.Error, Pvalue = `p-value`)]


alloutput1 <- dplyr::left_join(alloutput, maineffectTable, by = c("competitionModel", "Species", "Climate")) %>%
  data.table
alloutput1[is.na(mainEffect), ':='(mainEffect = 0, mainEffect_SE = 0)]

alloutput1 <- dplyr::left_join(alloutput1, interactionTable, by = c("competitionModel", "Species", "Climate")) %>%
  data.table
alloutput1 <- alloutput1[is.na(interactEff),':='(interactEff = 0, 
                                                 interactEff_SE = 0)]
alloutput1[,':='(EffectMin = mainEffect+CompetitionMinctd*interactEff, 
                 EffectMin_SE = sqrt(mainEffect_SE^2+interactEff_SE^2),
                 EffectMax = mainEffect+CompetitionMaxctd*interactEff, 
                 EffectMax_SE = sqrt(mainEffect_SE^2+interactEff_SE^2))]
alloutput1[, lineTransp:=1]
alloutput1[Pvalue>0.05,lineTransp:=2]
alloutput1[, mainPointShape:=1]
alloutput1[MainPvalue>0.05, mainPointShape:=2]

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
climateWithCompTable[,':='(Species = factor(Species, levels = c("All trees", majorspecies[2:4], 
                                                                "Minor species group")))]

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
ZeroLines[,':='(x = -Inf, xend = Inf, y = 0, yend = 0,
                SeasonComp = factor(SeasonComp, levels = c("Annual anomaly", 
                                                           "Growing season anomaly", 
                                                           "Non-growing season anomaly")),
                ClimateName = factor(ClimateName, levels = c("Temperature", "CMI", "CO2")))]
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
newlabels1 <- list("Temperature" = expression(atop("Sensitivity to temperature anomaly",
                                                   paste("(", cm^{2}, " ", year^{-1}, " ", degree, C^{-1}, " per tree)"))),
                   "CMI" = expression(atop("Sensitivity to CMI anomaly",
                                           paste("(", cm^{2}, " ", year^{-1}, " ", mm^{-1}, " per tree)"))),
                   'CO2'= expression(atop(paste("Sensitivity to C", O[2], " anomaly"),
                                          paste("(", cm^{2}, " ", year^{-1}, " ", ppm^{-1}, " per tree)"))))
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

climateWithCompTable[EffectMin>=EffectMax,
                     ':='(EffectStart = EffectMin+1.98*EffectMin_SE,
                          EffectEnd = EffectMax-1.98*EffectMax_SE)]
climateWithCompTable[EffectMin<EffectMax,
                     ':='(EffectStart = EffectMin-1.98*EffectMin_SE,
                          EffectEnd = EffectMax+1.98*EffectMax_SE)]
climateWithCompTable[,':='(EffectStart = exp(EffectStart)-1,
                           EffectEnd = exp(EffectEnd)-1)]
FigureB <- ggplot(data = climateWithCompTable[lineTransp == 1,], 
                  aes(x = xaxis, y = mainEffect))+
  facet_grid(ClimateName~SeasonComp,
             scales = "free_y", switch = "y",
             labeller = figure_labeller, drop = FALSE)+
  geom_segment(data = ZeroLines, aes(x = x, xend = xend, 
                                     y = y, yend = yend), linetype = 2, col = "gray", size = 1)+
  geom_segment(aes(col = Species, y = EffectStart, 
                   yend = EffectEnd, 
                   x = xaxis, xend = xaxis), 
               arrow = arrow(length = unit(0.05, "npc")), size = 1)+
  geom_errorbar(data = climateWithCompTable[lineTransp == 2,],
                aes(col = Species, ymin = pmin(EffectStart), 
                    ymax = pmax(EffectEnd), 
                    x = xaxis), 
                size = 1, width = 0.2)+
  scale_x_continuous(name = "a", limits = c(0.5, 5.5))+
  geom_segment(data = Yaxislines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1.5, col = "black")+
  geom_segment(data = Xaxislines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1, col = "black")+
  geom_text(data = labeltexts, aes(x = x, y = y, label = labels),
            vjust = 1.5, hjust = -1.5, size = 9)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85, 0.50),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.background = element_rect(colour = "black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 17, face = "italic"))
workPath <- "~/GitHub/Climate_Growth"

if(useLogQunatile95){
  ggsave(file.path(workPath, "TablesFigures",
                   "Figure S8. Climate associations with H_BA.png"), FigureB,
         width = 13, height = 11)
} else {
  ggsave(file.path(workPath, "TablesFigures",
                   "Figure 3. Climate associations with H.png"), FigureB,
         width = 10, height = 10)
}



