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
  if(indispecies == "All"){
    speciesData <- data.table::copy(analysesData)
  } else if(indispecies == "Other"){
    speciesData <- analysesData[!(Species %in% majorSpecies),]
  } else {
    speciesData <- analysesData[Species == indispecies,]
  }
  climateWithClimateTable <- rbind(data.table(expand.grid(Species = indispecies,
                                                              Climate = as.character(longcol), 
                                                              CompetitionName = "IntraH",
                                                              Competition = c(min(speciesData$IntraH1_3),
                                                                         exp(mean(log(speciesData$IntraH1_3+1)))-1,
                                                                         max(speciesData$IntraH1_3)),
                                                              stringsAsFactors = FALSE)),
                                       data.table(expand.grid(Species = indispecies,
                                                              Climate = as.character(longcol), 
                                                              CompetitionName = "InterH",
                                                              Competition = c(min(speciesData$InterH0_4),
                                                                         exp(mean(log(speciesData$InterH0_4+1)))-1,
                                                                         max(speciesData$InterH0_4)),
                                                              stringsAsFactors = FALSE)),
                                       data.table(expand.grid(Species = indispecies,
                                                              Climate = as.character(longcol), 
                                                              CompetitionName = "RBI",
                                                              Competition = c(min(speciesData$RBI),
                                                                      mean(speciesData$RBI),
                                                                      max(speciesData$RBI)),
                                                              stringsAsFactors = FALSE)),
                                       data.table(expand.grid(Species = indispecies,
                                                              Climate = as.character(longcol), 
                                                              CompetitionName = "SB",
                                                              Competition = c(min(speciesData$PlotBiomass),
                                                                     exp(mean(log(speciesData$PlotBiomass))),
                                                                     max(speciesData$PlotBiomass)),
                                                              stringsAsFactors = FALSE))) %>% data.table
  climateWithClimateTable[CompetitionName == "IntraH",
                          ':='(Competitionctd = log(Competition+1)-mean(log(speciesData$IntraH1_3+1)))]
  climateWithClimateTable[CompetitionName == "InterH",
                          ':='(Competitionctd = log(Competition+1)-mean(log(speciesData$InterH0_4+1)))]
  climateWithClimateTable[CompetitionName == "RBI",
                          ':='(Competitionctd = Competition-mean(speciesData$RBI))]
  climateWithClimateTable[CompetitionName == "SB",
                          ':='(Competitionctd = log(Competition)-mean(log(speciesData$PlotBiomass)))]
  climateWithClimateTable[,xscale:=0]
  
  climateWithClimateTable[CompetitionName == "IntraH" & Competition == min(speciesData$IntraH1_3), 
                          xscale:=as.numeric(m)]
  climateWithClimateTable[CompetitionName == "IntraH" & Competition == exp(mean(log(speciesData$IntraH1_3+1)))-1, 
                          xscale:=as.numeric(m+1)]
  climateWithClimateTable[CompetitionName == "IntraH" & Competition == max(speciesData$IntraH1_3),
                          xscale:=as.numeric(m+2)]
  climateWithClimateTable[CompetitionName == "InterH" & Competition == min(speciesData$InterH0_4), 
                          xscale:=as.numeric(m)]
  climateWithClimateTable[CompetitionName == "InterH" & Competition == exp(mean(log(speciesData$InterH0_4+1)))-1,
                          xscale:=as.numeric(m+1)]
  climateWithClimateTable[CompetitionName == "InterH" & Competition == max(speciesData$InterH0_4),
                          xscale:=as.numeric(m+2)]
  climateWithClimateTable[CompetitionName == "RBI" & Competition == min(speciesData$RBI),
                          xscale:=as.numeric(m+2)]
  climateWithClimateTable[CompetitionName == "RBI" & Competition == mean(speciesData$RBI),
                          xscale:=as.numeric(m+1)]
  climateWithClimateTable[CompetitionName == "RBI" & Competition == max(speciesData$RBI), 
                          xscale:=as.numeric(m)]
  climateWithClimateTable[CompetitionName == "SB" & Competition == min(speciesData$PlotBiomass), 
                          xscale:=as.numeric(m)]
  climateWithClimateTable[CompetitionName == "SB" & Competition == exp(mean(log(speciesData$PlotBiomass))), 
                          xscale:=as.numeric(m+1)]
  climateWithClimateTable[CompetitionName == "SB" & Competition == max(speciesData$PlotBiomass),
                          xscale:=as.numeric(m+2)]
  m <- m+4

  if(indispecies == "All"){
    alloutput <- climateWithClimateTable
  } else {
    alloutput <- rbind(alloutput, climateWithClimateTable)
  }
}


maineffectTable <- OverallResults[rn == "Climatectd",][, .(Species, Climate, mainEffect = Value,
                                                           mainEffect_SE = Std.Error)]

interactionTable <- OverallResults[rn %in% c("Climatectd:logIntraHctd", "Climatectd:logInterHctd",
                                             "Climatectd:RBIctd", "Climatectd:logSBctd" ),][
                                               ,.(rn, Species, Climate, interactEff = Value, 
                                                  interactEff_SE = Std.Error)]
interactionTable[rn == "Climatectd:logIntraHctd", CompetitionName:="IntraH"]
interactionTable[rn == "Climatectd:logInterHctd", CompetitionName:="InterH"]
interactionTable[rn == "Climatectd:RBIctd", CompetitionName:="RBI"]
interactionTable[rn == "Climatectd:logSBctd", CompetitionName:="SB"]
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
climateWithCompTable[Climate %in% c("ATA", "ACMIA", "ACO2A") & CompetitionName == "RBI",
                     SeasonComp:="WholeRBI"]
climateWithCompTable[Climate %in% c("ATA", "ACMIA", "ACO2A") & CompetitionName == "SB",
                     SeasonComp:="WholeSB"]

climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A") & CompetitionName == "IntraH",
                     SeasonComp:="GSIntraH"]
climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A") & CompetitionName == "InterH",
                     SeasonComp:="GSInterH"]
climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A") & CompetitionName == "RBI",
                     SeasonComp:="GSRBI"]
climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A") & CompetitionName == "SB",
                     SeasonComp:="GSSB"]

climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A") & CompetitionName == "IntraH", 
                     SeasonComp:="NGSIntraH"]
climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A") & CompetitionName == "InterH", 
                     SeasonComp:="NGSInterH"]
climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A") & CompetitionName == "RBI", 
                     SeasonComp:="NGSRBI"]
climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A") & CompetitionName == "SB", 
                     SeasonComp:="NGSSB"]

climateWithCompTable[, SeasonComp:=factor(SeasonComp, 
                                          levels = c("WholeIntraH", "WholeInterH", "WholeRBI", "WholeSB",
                                                     "GSIntraH", "GSInterH", "GSRBI", "GSSB",
                                                     "NGSIntraH", "NGSInterH", "NGSRBI", "NGSSB"))]

climateWithCompTable[,':='(Species = factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                                           labels = c("All species", "Jack pine", "Trembling aspen",
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

mainEffect <- climateWithCompTable[xscale %in% c(2, 6, 10, 14, 18) & mainEffect != 0,][,':='(x = xscale)]


newlabels1 <- list("Temperature" = "Effect of temperature",
                   "CMI" = "Effect of CMI",
                   'CO2'= expression(paste("Effect of C", O[2])))
newlabels2 <- list("WholeIntraH" = "IntraH",
                   "WholeInterH" = "InterH",
                   "WholeRBI" = "RBI",
                   "WholeSB" = "SB",
                   "GSIntraH" = "IntraH",
                   "GSInterH" = "InterH",
                   "GSRBI" = "RBI",
                   "GSSB" = "SB",
                   "NGSIntraH" = "IntraH",
                   "NGSInterH" = "InterH",
                   "NGSRBI" = "RBI",
                   "NGSSB" = "SB")


figure_labeller <- function(variable,value){
  if(variable == "ClimateName"){
     return(newlabels1[value])
  } else if (variable == "SeasonComp"){
    return(newlabels2[value])
  } 
}
majorYpanelbreaklines <- data.table(SeasonComp = factor(c("GSIntraH", "NGSIntraH"),
                                                        levels = c("WholeIntraH", "WholeInterH", "WholeRBI", "WholeSB",
                                                                   "GSIntraH", "GSInterH", "GSRBI", "GSSB",
                                                                   "NGSIntraH", "NGSInterH", "NGSRBI", "NGSSB")),
                                    x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

minorYpanelbreaklines <- data.table(SeasonComp = factor(c("WholeInterH", "WholeRBI", "WholeSB",
                                                          "GSInterH", "GSRBI", "GSSB",
                                                          "NGSInterH", "NGSRBI", "NGSSB"),
                                                        levels = c("WholeIntraH", "WholeInterH", "WholeRBI", "WholeSB",
                                                                   "GSIntraH", "GSInterH", "GSRBI", "GSSB",
                                                                   "NGSIntraH", "NGSInterH", "NGSRBI", "NGSSB")),
                                    x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

figureB <- ggplot(data = segmentPoints[Effect != Effectend, ], aes(x = x, y = Effect))+
  geom_segment(aes(x = -1, xend = 20, y = 0, yend = 0), linetype = 2, col = "gray", size = 1)+
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
               size = 3, col = "white")+
  geom_segment(data = minorYpanelbreaklines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1, col = "white")+
  geom_point(data = mainEffect[interactEff == 0, ], aes(x = x, y = mainEffect, group = Species, col = Species))+
  facet_grid(ClimateName~SeasonComp,
             scales = "free_y", switch = "both", 
             labeller = figure_labeller, drop = FALSE)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.margin.x = unit(0, "lines"), 
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        legend.position = c(0.8, 0.6),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.background = element_rect(colour = "black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 15, vjust = -2))





climates <- read.csv(file.path(workPath, "data", "plotClimates.csv"), header = TRUE,
                          stringsAsFactors = FALSE) %>%
  data.table
climates[,Year:=(FinYear+IniYear)/2]
climate_longform <- reshape(data = climates, varying = longcol, v.names = "Value",
                            times = longcol, timevar = "DependentVariable", 
                            direction = "long") %>% data.table
climate_longform[,Yearctd:=Year-mean(Year)]


for(indiclimategroup in c("temperature", "CMI", "CO2")){
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
allClimateData[, Season:=factor(Season, levels = c("Whole year", "Growing season", "Non-growing season"))]

textYdif <- 0.3
textYcentral <- -1.1
Figure_a <- ggplot(data = allClimateData[Climate %in% temperature],
                            aes(x = Year, y = Value))+
  geom_point(aes(col = Season), alpha = 0.1)+
  geom_line(aes(x = Year, y = predValue, col = Season, linetype = linetype), size = 1)+
  scale_linetype(guide = "none")+
  scale_color_manual(name = "Season", 
                     values = c("red", "green", "blue"),
                     label = c("Whole year", "Growing season", "Non-growing season"))+
  scale_y_continuous(name = expression(paste("Temperature anomaly (", degree, "C)")), limits = c(-1.5, 1.5),
                     breaks = seq(-1.5, 1.5, by = 0.5))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  annotate("text", x = c(1985, 1999, rep(2002, 3)), 
           y = c(1.5, textYcentral, textYcentral+textYdif, textYcentral, textYcentral-textYdif),
           label = c("a", "Slope:", as.character(allcoeff[DependentVariable %in% temperature,]$Value)),
           size = c(10, rep(5, 4)), col = c("black", "black", "red", "green", "blue"))+
  annotate("text", x = 2007, y = textYcentral, label = paste("~degree~C~year^{-1}"),
           size = 5, parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = c(0.3, 0.82),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))



# the a
Figure_b <- ggplot(data = climateWithRBITable[Climate %in% temperature,], 
                aes(x = RBI, y = climateEffect))+
  geom_segment(x = 0, xend = 100, y = 0, yend = 0, 
               linetype = 2, size = 1, colour = "gray")+
  geom_line(aes(group = interaction(Species, Climate), col = Species), size = 1)+
  geom_point(data = maineffectTable[Climate %in% temperature,],
             aes(x = RBI, y = climateEffect, col = Species), size = 2)+
  geom_errorbar(data = maineffectTable[Climate %in% temperature,],
                aes(x = RBI, ymin = climateEffect-1.98*maineffect_SE, 
                    ymax = climateEffect+1.98*maineffect_SE, col = Species))+
  facet_wrap(~Season)+
  geom_rect(aes(xmin = 15, xmax = 85, ymin = 1, ymax = 1.2), col = "gray", fill = "gray")+
  geom_text(data = data.frame(RBI = 0, climateEffect = 1.2, label = "b", Season = "Whole year"),
            aes(x = RBI, y = climateEffect, label = label), size = 10)+
  geom_text(data = data.frame(RBI = 50, climateEffect = 1.1, 
                                label = c("Whole year", "Growing season", "Non-growing season"),
                               Season = c("Whole year", "Growing season", "Non-growing season")),
              aes(x = RBI, y = climateEffect, label = label), size = 5)+
  scale_color_manual(name = "Species", 
                     values = c("black", "red", "green", "blue", "yellow"),
                     labels = c("All species", "Jack pine", "Trembling aspen",
                                "Black spruce", "Other species"))+
  scale_linetype(guide = "none")+
  scale_y_continuous(name = "Temperature effect",
                     limits = c(-0.4, 1.2), 
                     breaks = round(seq(-0.4, 1.2, by = 0.4),1))+
  scale_x_continuous(name = "RBI", limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # axis.line.x = element_line(colour = "black", size = 1),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black", size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = "none")



PrepFigure <- "FALSE"
if(PrepFigure){
  print(allcoeff[DependentVariable %in% precipitationNames,.(DependentVariable, Value)])
  precipitationFigure_c <- ggplot(data = allClimateData[DependentVariable %in% precipitationNames],
                                  aes(x = Year, y = Value))+
    geom_point(aes(group = interaction(PlotID, DependentVariable), col = DependentVariable),
               alpha = 0.1)+
    geom_line(aes(x = Year, y = predValue, col = DependentVariable, linetype = linetype), size = 1)+
    scale_linetype(guide = "none")+
    scale_color_manual(name = "Season", 
                       values = c("red", "green", "blue"),
                       label = c(expression(paste("Whole year (slope: 3.65 mm ", year^{-1},")")),
                                 expression(paste("Growing season (slope: 2.86 mm ", year^{-1},")")),
                                 expression(paste("Non-growing season (slope: 1.04 mm ", year^{-1},")"))))+
    scale_y_continuous(name = paste("Precipitation anomaly (mm)"), limits = c(-100, 130),
                       breaks = seq(-100, 130, by = 35))+
    scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
    annotate("text", x = 1985, y = 130, label = "c", size = 10)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_blank(),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15))
  
  
  precipitationFigure_d <- ggplot(data = climateWithRBITable[Climate %in% precipitationNames,], 
                                  aes(x = newDomincneIndex, y = climateEffect))+
    geom_segment(data = data.frame(x = c(0, 110, 220), xend = c(100, 210, 320)),
                 aes(x = x, xend = xend, y = 0, yend = 0), 
                 linetype = 2, size = 1, colour = "gray")+
    geom_line(aes(group = interaction(Species, Climate), col = Species, linetype = linetype), size = 1)+
    geom_point(data = maineffectTable[Climate %in% precipitationNames,],
               aes(x = newDomincneIndex, y = climateEffect, col = Species), size = 2)+
    geom_errorbar(data = maineffectTable[Climate %in% precipitationNames,],
                  aes(x = newDomincneIndex, ymin = climateEffect-1.98*maineffect_SE, 
                      ymax = climateEffect+1.98*maineffect_SE, col = Species))+
    scale_color_manual(name = "Species", 
                       values = c("red", "green", "blue"),
                       labels = c("Jack pine", "Trembling aspen", "Black spruce"))+
    scale_linetype(guide = "none")+
    scale_y_continuous(name = "Effect of precipitation on growth",
                       limits = c(-0.025, 0.015), 
                       breaks = round(seq(-0.02, 0.01, by = 0.1),2))+
    scale_x_continuous(name = "RBI index", limits = c(0, 320),
                       breaks = c(seq(0, 100, by = 20), seq(110, 210, by = 20),
                                  seq(220, 320, 20)))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15))
}


print(allcoeff[DependentVariable %in% CMI,.(DependentVariable, Value)])
textYdif <- 30
textYcentral <- -100
Figure_c <- ggplot(data = allClimateData[Climate %in% CMI,],
                              aes(x = Year, y = Value))+
  geom_point(aes(col = Season),
             alpha = 0.1)+
  geom_line(aes(x = Year, y = predValue, col = Season, linetype = linetype), size = 1)+
  scale_linetype(guide = "none")+
  scale_color_manual(name = "Season", 
                     values = c("red", "green", "blue"),
                     label = c(expression(paste("Whole year (slope: 4.18 mm ", year^{-1},")")),
                               expression(paste("Growing season (slope: 4.18 mm ", year^{-1},")")),
                               expression(paste("Non-growing season (slope: 0.31 mm ", year^{-1},")"))))+
  scale_y_continuous(name = paste("CMI anomaly (mm)"),
                     limits = c(-140, 130),
                     breaks = seq(-140, 130, by = 40))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  annotate("text", x = c(1999, rep(2002, 3)), 
           y = c(textYcentral, textYcentral+textYdif, textYcentral, textYcentral-textYdif),
           label = c("Slope:", as.character(allcoeff[DependentVariable %in% CMI,]$Value)),
           size = c(rep(5, 4)), col = c("black", "red", "green", "blue"))+
  annotate("text", x = 2007, y = textYcentral, label = paste("~mm~year^{-1}"),
           size = 5, parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))


CMINames1 <- CMI[1:2]
climateWithRBITable[Climate %in% CMINames1, climateEffect:=climateEffect*100]
maineffectTable[Climate %in% CMINames1,':='(climateEffect = climateEffect*100,
                                            maineffect_SE = maineffect_SE*100)]
figuredatacmi <- climateWithRBITable[Climate %in% CMI,][
  !(Season == "Non-growing season" & RBI > 60), ][
    !(Season == "Non-growing season" & RBI < 40),]
  
Figure_d <- ggplot(data = figuredatacmi, aes(x = RBI, y = climateEffect))+
  geom_segment(data = data.frame(Climate = CMINames1, x = 0, xend = 100, y = 0, yend = 0),
               aes(x = x, xend = xend, y = y, yend = yend), 
               linetype = 2, size = 1, colour = "gray")+
  geom_line(aes(col = Species), size = 1)+
  facet_wrap(~Climate, nrow = 1)+
  geom_point(data = maineffectTable[Climate %in% CMINames1,],
             aes(x = RBI, y = climateEffect, col = Species), size = 2)+
  geom_errorbar(data = maineffectTable[Climate %in% CMINames1,],
                aes(x = RBI, ymin = climateEffect-1.98*maineffect_SE, 
                    ymax = climateEffect+1.98*maineffect_SE, col = Species))+
  scale_color_manual(name = "Species", 
                      values = c("black", "red", "green", "blue", "yellow"),
                      labels = c("All species", "Jack pine", "Trembling aspen",
                                 "Black spruce", "Other species"))+
  scale_y_continuous(name = expression(paste("CMI effect (",10^{-2}, ")")),
                     limits = c(-1.03, 0.8), 
                     breaks = round(seq(-1, 0.8, by = 0.3),1))+
  scale_x_continuous(name = "RBI index", limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = c(0.8, 0.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))


print(allcoeff[DependentVariable %in% CO2,.(DependentVariable, Value)])
textYdif <- 4
textYcentral <- -12
Figure_e <- ggplot(data = allClimateData[Climate %in% CO2,],
                    aes(x = Year, y = Value))+
  geom_point(aes(col = Season),
             alpha = 0.1)+
  geom_line(aes(x = Year, y = predValue, col = Season, linetype = linetype), size = 1)+
  scale_linetype(guide = "none")+
  scale_color_manual(name = "Season", 
                     values = c("red", "green", "blue"),
                     label = c(expression(paste("Whole year (slope: 1.80 ppm ", year^{-1},")")),
                               expression(paste("Growing season (slope: 1.80 ppm ", year^{-1},")")),
                               expression(paste("Non-growing season (slope: 1.79 ppm ", year^{-1},")"))))+
  scale_y_continuous(name = expression(paste(CO[2], " anomaly (ppm)")), limits = c(-20, 22),
                     breaks = seq(-20, 20, by = 8))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  annotate("text", x = c(1999, rep(2002, 3)), 
           y = c(textYcentral, textYcentral+textYdif, textYcentral, textYcentral-textYdif),
           label = c("Slope:", "1.80", "1.80", "1.79"),
           size = c(rep(5, 4)), col = c("black", "red", "green", "blue"))+
  annotate("text", x = 2007, y = textYcentral, label = paste("~ppm~year^{-1}"),
           size = 5, parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "none",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))
climateWithRBITable[Climate %in% CO2,climateEffect:=climateEffect*100]
maineffectTable[Climate %in% CO2,':='(climateEffect = climateEffect*100, maineffect_SE = maineffect_SE*100)]

Figure_f <- ggplot(data = climateWithRBITable[Climate %in% CO2,], 
                      aes(x = RBI, y = climateEffect))+
  geom_segment(x = 0, xend = 100, y = 0, yend = 0, 
               linetype = 2, size = 1, colour = "gray")+
  geom_line(aes(col = Species), size = 1)+
  facet_wrap(~Climate, nrow = 1)+
  geom_point(data = maineffectTable[Climate %in% CO2,],
             aes(x = RBI, y = climateEffect, col = Species), size = 2)+
  geom_errorbar(data = maineffectTable[Climate %in% CO2,],
                aes(x = RBI, ymin = climateEffect-1.98*maineffect_SE, 
                    ymax = climateEffect+1.98*maineffect_SE, col = Species))+
  scale_color_manual(name = "Species", 
                     values = c("black", "red", "green", "blue", "yellow"),
                     labels = c("All species", "Jack pine", "Trembling aspen",
                                "Black spruce", "Other species"))+
  scale_y_continuous(name = expression(paste(CO[2], " effect (", 10^{-2}, ")")),
                     limits = c(-5.5, 2), 
                     breaks = round(seq(-5.5, 2, by = 1.5),2))+
  scale_x_continuous(name = "RBI", limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "none")



Figa_Grob <- ggplotGrob(Figure_a)
Figb_Grob <- ggplotGrob(Figure_b)
Figc_Grob <- ggplotGrob(Figure_c)
Figd_Grob <- ggplotGrob(Figure_d)
Fige_Grob <- ggplotGrob(Figure_e)
Figf_Grob <- ggplotGrob(Figure_f)

# height
Figa_Grob$heights <- Fige_Grob$heights
Figb_Grob$heights <- Figf_Grob$heights
Figc_Grob$heights <- Fige_Grob$heights
Figd_Grob$heights <- Figf_Grob$heights

# width
Figa_Grob$widths <- Figc_Grob$widths
Fige_Grob$widths <- Figc_Grob$widths


Figb_Grob$widths <- Figf_Grob$widths
Figd_Grob$widths <- Figf_Grob$widths


dev(4)
clearPlot()
plotlayout <- rbind(c(1, 2, 2), c(3, 4, 4), c(5, 6, 6))
c <- grid.arrange(Figa_Grob, Figb_Grob, Figc_Grob,Figd_Grob, Fige_Grob, Figf_Grob,
                  layout_matrix = plotlayout)

ggsave(file = file.path(workPath, "TablesFigures", "ClimateTrandsAndEffects.png"), c,
       width = 18, height = 10)


