rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "ClimateModelSelection.RData"))
rm(speciesData, indispecies, fullthemodel, i, indiclimate, newSigNIDV, reducedFormu,
   reducedModel, signIDV)

for(i in 1:length(allFixedCoeff)){
  modelName <- unlist(strsplit(names(allFixedCoeff)[i], "_", fixed = T))
  coeffIndi <- cbind(data.table(Model = rep(names(allFixedCoeff)[i], nrow(allFixedCoeff[[i]])),
                                Species = modelName[1],
                                Climate = modelName[2]), 
                     
                     allFixedCoeff[[i]])
  if(i == 1){
    allCoeff <- coeffIndi
  } else {
    allCoeff <- rbind(allCoeff, coeffIndi)
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

m <- 1
for(indispecies in studySpecies){
  
  speciesData <- analysesData[DataType == indispecies,]
  
  climateWithClimateTable <- rbind(data.table(expand.grid(Species = indispecies,
                                                          Climate = as.character(longcol), 
                                                          CompetitionName = "IntraH",
                                                          Competition = c(min(speciesData$IntraH),
                                                                          exp(mean(log(speciesData$IntraH+1)))-1,
                                                                          max(speciesData$IntraH)),
                                                          stringsAsFactors = FALSE)),
                                   data.table(expand.grid(Species = indispecies,
                                                          Climate = as.character(longcol), 
                                                          CompetitionName = "InterH",
                                                          Competition = c(min(speciesData$InterH),
                                                                          exp(mean(log(speciesData$InterH+1)))-1,
                                                                          max(speciesData$InterH)),
                                                          stringsAsFactors = FALSE))) %>% data.table
  climateWithClimateTable[CompetitionName == "IntraH",
                          ':='(Competitionctd = log(Competition+1)-mean(log(speciesData$IntraH+1)))]
  climateWithClimateTable[CompetitionName == "InterH",
                          ':='(Competitionctd = log(Competition+1)-mean(log(speciesData$InterH+1)))]
  climateWithClimateTable[,xscale:=0]
  
  climateWithClimateTable[CompetitionName == "IntraH" & Competition == min(speciesData$IntraH), 
                          xscale:=as.numeric(m)]
  climateWithClimateTable[CompetitionName == "IntraH" & Competition == exp(mean(log(speciesData$IntraH+1)))-1, 
                          xscale:=as.numeric(m+1)]
  climateWithClimateTable[CompetitionName == "IntraH" & Competition == max(speciesData$IntraH),
                          xscale:=as.numeric(m+2)]
  climateWithClimateTable[CompetitionName == "InterH" & Competition == min(speciesData$InterH), 
                          xscale:=as.numeric(m)]
  climateWithClimateTable[CompetitionName == "InterH" & Competition == exp(mean(log(speciesData$InterH+1)))-1,
                          xscale:=as.numeric(m+1)]
  climateWithClimateTable[CompetitionName == "InterH" & Competition == max(speciesData$InterH),
                          xscale:=as.numeric(m+2)]
  m <- m+4
  
  if(indispecies == "All species"){
    alloutput <- climateWithClimateTable
  } else {
    alloutput <- rbind(alloutput, climateWithClimateTable)
  }
}


maineffectTable <- OverallResults[rn == "Climatectd",][, .(Species, Climate, mainEffect = Value,
                                                           mainEffect_SE = Std.Error)]

interactionTable <- OverallResults[rn %in% c("Climatectd:logIntraHctd", "Climatectd:logInterHctd"),][
  ,.(rn, Species, Climate, interactEff = Value, 
     interactEff_SE = Std.Error)]
interactionTable[rn == "Climatectd:logIntraHctd", CompetitionName:="IntraH"]
interactionTable[rn == "Climatectd:logInterHctd", CompetitionName:="InterH"]
interactionTable[,rn:=NULL]

alloutput1 <- dplyr::left_join(alloutput, maineffectTable, by = c("Species", "Climate")) %>%
  data.table
alloutput1[is.na(mainEffect), ':='(mainEffect = 0, mainEffect_SE = 0)]

alloutput1 <- dplyr::left_join(alloutput1, interactionTable, by = c("Species", "Climate", "CompetitionName")) %>%
  data.table
alloutput1 <- alloutput1[is.na(interactEff),':='(interactEff = 0, 
                                                 interactEff_SE = 0)]
alloutput1[,':='(Effect = mainEffect+Competitionctd*interactEff, 
                 Effect_SE = sqrt(mainEffect_SE^2+interactEff_SE^2))]
alloutput1[xscale %in% c(2, 6, 10, 14, 18), ':='(Effect = mainEffect, 
                                                 Effect_SE = mainEffect_SE)]
climateWithCompTable <- data.table::copy(alloutput1)
climateWithCompTable[Climate %in% c("ATA", "ACMIA", "ACO2A") & CompetitionName == "IntraH",
                     SeasonComp:="WholeIntraH"]
climateWithCompTable[Climate %in% c("ATA", "ACMIA", "ACO2A") & CompetitionName == "InterH",
                     SeasonComp:="WholeInterH"]

climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A") & CompetitionName == "IntraH",
                     SeasonComp:="GSIntraH"]
climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A") & CompetitionName == "InterH",
                     SeasonComp:="GSInterH"]

climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A") & CompetitionName == "IntraH", 
                     SeasonComp:="NGSIntraH"]
climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A") & CompetitionName == "InterH", 
                     SeasonComp:="NGSInterH"]

climateWithCompTable[, SeasonComp:=factor(SeasonComp, 
                                          levels = c("WholeIntraH", "WholeInterH", 
                                                     "GSIntraH", "GSInterH", 
                                                     "NGSIntraH", "NGSInterH"))]

climateWithCompTable[,':='(Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen",
                                                                "Black spruce", "Other species")))]

climateWithCompTable[Climate %in% c("ATA", "GSTA", "NONGSTA"), ClimateName:="Temperature"]
climateWithCompTable[Climate %in% c("ACMIA", "GSCMIA", "NONGSCMIA"), ClimateName:="CMI"]
climateWithCompTable[Climate %in% c("ACO2A", "GSCO2A", "NONGSCO2A"), ClimateName:="CO2"]
climateWithCompTable[, ClimateName:=factor(ClimateName, levels = c("Temperature", "CMI", "CO2"))]



startPoints <- climateWithCompTable[xscale %in% c(1, 5, 9, 13, 17),][, ':='(x = xscale+1)]
endPoints <- climateWithCompTable[xscale %in% c(3, 7, 11, 15, 19), ][, .(Species, Climate, CompetitionName, 
                                                                         xend = xscale-1, 
                                                                         Effectend = Effect, Effectend_SE = Effect_SE)]
segmentPoints <- setkey(startPoints, Species, Climate, CompetitionName)[setkey(endPoints, Species, Climate, CompetitionName), 
                                                                        nomatch = 0]
segmentPoints[,':='(Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen",
                                                         "Black spruce", "Other species")))]
mainEffect <- climateWithCompTable[xscale %in% c(2, 6, 10, 14, 18) & mainEffect != 0,][,':='(x = xscale)]
mainEffect[,':='(Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen",
                                                         "Black spruce", "Other species")))]

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
majorYpanelbreaklines <- data.table(SeasonComp = factor(c("WholeIntraH", "GSIntraH", "NGSIntraH"),
                                                        levels = c("WholeIntraH", "WholeInterH", 
                                                                   "GSIntraH", "GSInterH", 
                                                                   "NGSIntraH", "NGSInterH")),
                                    x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

majorXpanelbreaklines <- data.table(ClimateName = factor(c("CO2"),
                                                        levels = c("Temperature", "CMI", 
                                                                   "CO2")),
                                    x = Inf, xend = -Inf, y = -Inf, yend = -Inf)
a <- unique(segmentPoints[,.(SeasonComp, ClimateName)],
            by = c("SeasonComp", "ClimateName"))
a[,':='(x = -1, xend = 20, Effect = 0, yend = 0)]
labeltexts <- data.table(ClimateName = factor("Temperature",
                                              levels = c("Temperature", "CMI", 
                                                         "CO2")), 
                         SeasonComp = factor(c("WholeIntraH", 
                                        "GSIntraH", "NGSIntraH"),
                                        levels = c("WholeIntraH", "WholeInterH", 
                                                   "GSIntraH", "GSInterH", 
                                                   "NGSIntraH", "NGSInterH")),
                         labels = letters[1:3],
                         x = 0.5, y = Inf)

FigureB <- ggplot(data = segmentPoints[Effect != Effectend, ], aes(x = x, y = Effect))+
  facet_grid(ClimateName~SeasonComp,
             scales = "free_y", switch = "both", 
             labeller = figure_labeller, drop = FALSE)+
  geom_point(data = mainEffect, aes(x = x, y = mainEffect, group = Species, col = Species),
             size = 0.5)+
  geom_segment(data = a, aes(x = x, xend = xend, 
                             y = Effect, yend = yend), linetype = 2, col = "gray", size = 1)+
  geom_errorbar(aes(group = Species, ymin = Effect-1.98*Effect_SE, ymax = Effect+1.98*Effect_SE), 
                width = 2, col = "gray", size = 1)+
  geom_errorbar(aes(group = Species, ymin = Effectend-1.98*Effectend_SE, ymax = Effectend+1.98*Effectend_SE), 
                width = 2, col = "gray", size = 1)+
  geom_errorbar(data = mainEffect[interactEff == 0, ], aes(x = x, ymin = mainEffect-mainEffect_SE, 
                                                           ymax = mainEffect+mainEffect_SE, group = Species),
                col = "gray", width = 2, size = 1)+
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
                 "Figure 5. Climate associations with tree growth.png"), FigureB,
       width = 10, height = 7)


temperature <- c("ATA", "GSTA", "NONGSTA")
allCMI <- c("ACMIA", "GSCMIA", "NONGSCMIA")
allCO2 <- c("ACO2A", "GSCO2A", "NONGSCO2A")
newlongcol <- c(temperature, allCMI, allCO2)
climates <- read.csv(file.path(workPath, "data", "plotClimates.csv"), header = TRUE,
                     stringsAsFactors = FALSE) %>%
  data.table
climates[,Year:=(FinYear+IniYear)/2]
climate_longform <- reshape(data = climates, varying = newlongcol, v.names = "Value",
                            times = newlongcol, timevar = "DependentVariable", 
                            direction = "long") %>% data.table
climate_longform[,Yearctd:=Year-mean(Year)]


for(indiclimategroup in c("temperature", "allCMI", "allCO2")){
  climateData <- climate_longform[DependentVariable %in% get(indiclimategroup),]
  climateModel <- lme(Value ~ DependentVariable/Yearctd, random =~(DependentVariable-1)|PlotID,
                      data = climateData)
  climateData$predValue <- predict(climateModel, newdata = climateData, level = 0) 
  coeff <- data.frame(summary(climateModel)$tTable)
  coeff <- coeff[row.names(coeff) %in% paste("DependentVariable", get(indiclimategroup), ":Yearctd", sep = ""),c(1,5)] %>%
    data.table
  coeff[,':='(DependentVariable = get(indiclimategroup),
              Value = round(coeff$Value, 2), linetype=1)]
  coeff[p.value>=0.05, linetype:=2]
  climateData <- setkey(climateData, DependentVariable)[setkey(coeff[,.(DependentVariable, linetype)], DependentVariable),
                                                        nomatch = 0]
  climateData <- climateData[,.(ClimateName = indiclimategroup, Climate = DependentVariable,
                                Year, Value, predValue, linetype)]
  if(indiclimategroup == "temperature"){ 
    allClimateData <- climateData
    allcoeff <- coeff
  } else {
    allClimateData <- rbind(allClimateData, climateData)
    allcoeff <- rbind(allcoeff, coeff)
  }
}
allClimateData$linetype <- factor(allClimateData$linetype, levels = c(1, 2))
allClimateData[Climate %in% c("ATA", "ACMIA", "ACO2A"), Season:="Whole year"]
allClimateData[Climate %in% c("GSTA", "GSCMIA", "GSCO2A"), Season:="Growing season"]
allClimateData[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A"), Season:="Non-growing season"]
allClimateData[, ':='(Season=factor(Season, levels = c("Whole year", "Growing season", "Non-growing season")),
                      ClimateName = factor(ClimateName, levels = c("temperature", "allCMI", "allCO2")))]
relativePosition <- 0.1
positionDiff <- 0.1

slopestexts1 <- allClimateData[,.(maxValue = max(Value), minValue = min(Value)), by = ClimateName]
slopestexts1[, ':='(y = minValue+(maxValue-minValue)*relativePosition, Year = 1999, labels = "Slope:")]
slopestexts2 <- allcoeff[,.(labels = as.character(Value), Climate = DependentVariable)]
slopestexts2[Climate %in% c("ATA", "ACMIA", "ACO2A"), Season:="Whole year"]
slopestexts2[Climate %in% c("GSTA", "GSCMIA", "GSCO2A"), Season:="Growing season"]
slopestexts2[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A"), Season:="Non-growing season"]
slopestexts2[Climate %in% temperature, ClimateName:="temperature"]
slopestexts2[Climate %in% allCMI, ClimateName:="allCMI"]
slopestexts2[Climate %in% allCO2, ClimateName:="allCO2"]
slopestexts2[, ':='(Season=factor(Season, levels = c("Whole year", "Growing season", "Non-growing season")),
                    ClimateName = factor(ClimateName, levels = c("temperature", "allCMI", "allCO2")))]
slopestexts2 <- dplyr::left_join(slopestexts2, slopestexts1[,.(ClimateName, maxValue, minValue)], by = "ClimateName") %>%
  data.table
slopestexts2[Season == "Whole year", y:=minValue+(maxValue-minValue)*(relativePosition+positionDiff)]
slopestexts2[Season == "Growing season", y:=minValue+(maxValue-minValue)*(relativePosition)]
slopestexts2[Season == "Non-growing season", y:=minValue+(maxValue-minValue)*(relativePosition-positionDiff)]
slopestexts2[, Year := 2002]
slopestexts3 <- data.table::copy(slopestexts1)[ClimateName == "temperature", UnitLabel := paste("~degree~C~year^{-1}")]
slopestexts3[ClimateName == "allCMI", UnitLabel := paste("~mm~year^{-1}")]
slopestexts3[ClimateName == "allCO2", UnitLabel := paste("~ppm~year^{-1}")]
slopestexts3[, Year:=2007]
slopestexts4 <- slopestexts3[ClimateName == "temperature", .(ClimateName, Year = 1985, y = Inf, label = "a")]

anewlaberller <- list("temperature" = expression(paste("Temperature anomaly (", degree, "C)")),
                   "allCMI" = "CMI anomaly (mm)",
                   'allCO2'= expression(paste("C", O[2], " anomaly (ppm)")))
figurea_labeller <- function(variable,value){
  if(variable == "ClimateName"){
    return(anewlaberller[value])
  } 
}

FigureA <- ggplot(data = allClimateData,
                   aes(x = Year, y = Value))+
  facet_grid(ClimateName~., switch = "y", scales = "free_y", labeller = figurea_labeller)+
  geom_point(aes(col = Season), alpha = 0.1)+
  geom_line(aes(x = Year, y = predValue, col = Season, linetype = linetype), size = 1)+
  scale_linetype(guide = "none")+
  scale_color_manual(name = "Season", 
                     values = c("red", "green", "blue"),
                     label = c("Whole year", "Growing season", "Non-growing season"))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  geom_text(data = slopestexts1, aes(x = Year, y = y, label = labels), size = 5)+
  geom_text(data = slopestexts2, aes(x = Year, y = y, label = labels, col = Season), 
            size = 5, show.legend = FALSE)+
  geom_text(data = slopestexts3, aes(x = Year, y = y, label = UnitLabel), size = 5, parse = TRUE)+
  # geom_text(data = slopestexts4, aes(x = Year, y = y, label = label), vjust = 1.5, size = 7)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 15, vjust = 2),
        strip.text = element_text(size = 15),
        strip.background = element_blank(),
        legend.title = element_text(size = 15),
        legend.position = c(0.2, 0.2),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size = 12))
ggsave(file.path(workPath, "TablesFigures", "Figure 4. regional climate changes.png"), FigureA,
       width = 7, height = 9)




