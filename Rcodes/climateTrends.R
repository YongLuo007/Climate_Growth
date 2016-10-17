rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "ClimateResults.RData"))
rm(speciesData, indispecies, fullthemodel, i, indiclimate, newSigNIDV, reducedFormu,
   reducedModel, signIDV)

tempallCoeff <- lapply(allClimateModels, function(x){data.table(summary(x)$tTable, keep.rownames = TRUE)})
for(i in 1:length(tempallCoeff)){
  modelName <- unlist(strsplit(names(tempallCoeff)[i], "_", fixed = T))
  coeffIndi <- cbind(data.table(Model = rep(names(tempallCoeff)[i], nrow(tempallCoeff[[i]])),
                                Species = modelName[1],
                                Climate = modelName[2]), 
                     
                     tempallCoeff[[i]])
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
CMI <- c("ACMIA", "GSCMIA", "NONGSCMIA")
CO2 <- c("ACO2A", "GSCO2A", "NONGSCO2A")
longcol <- c(temperature, CMI, CO2)
climateWithRBITable <- data.table(expand.grid(Species = as.character(studySpecies), 
                                              Climate = as.character(longcol), 
                                              RBI = seq(0, 100, by = 10),
                                              stringsAsFactors = FALSE))

maineffectTable <- data.table(Species = OverallResults[rn == "climatectd",]$Species,
                              Climate = OverallResults[rn == "climatectd",]$Climate,
                              maineffect = OverallResults[rn == "climatectd",]$Value,
                              maineffect_SE = OverallResults[rn == "climatectd",]$Std.Error)
interactionTable <- data.table(Species = OverallResults[rn == "climatectd:RBIctd",]$Species,
                               Climate = OverallResults[rn == "climatectd:RBIctd",]$Climate,
                               interactioneffect = OverallResults[rn == "climatectd:RBIctd",]$Value,
                               Pvalue = OverallResults[rn == "climatectd:RBIctd",]$'p-value')
climateWithRBITable <- dplyr::left_join(climateWithRBITable, maineffectTable, by = c("Species", "Climate")) %>%
  data.table
climateWithRBITable[is.na(maineffect), ':='(maineffect = 0, maineffect_SE = 0)]
climateWithRBITable <- dplyr::left_join(climateWithRBITable, interactionTable, by = c("Species", "Climate")) %>%
  data.table
climateWithRBITable <- climateWithRBITable[!is.na(interactioneffect),]


speciesdata <- rbind(copy(analysesData)[,.(RBI, Species="All")],
                     analysesData[,.(RBI, Species)][!(Species %in% majorSpecies), Species:="Other"])

speciesdata <- speciesdata[,.(meanRBI = mean(RBI)), by = Species]
climateWithRBITable <- setkey(climateWithRBITable, Species)[setkey(speciesdata, Species),
                                                                        nomatch = 0]
climateWithRBITable[,':='(climateEffect = (RBI-meanRBI)*interactioneffect+maineffect,
                          Species = factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                                           labels = c("All species", "Jack pine", "Trembling aspen",
                                                      "Black spruce", "Other species")))]
climateWithRBITable[Climate %in% c("ATA", "ACMIA", "ACO2A"), Season:="Whole year"]
climateWithRBITable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A"), Season:="Growing season"]
climateWithRBITable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A"), Season:="Non-growing season"]
climateWithRBITable[, Season:=factor(Season, levels = c("Whole year", "Growing season", "Non-growing season"))]
climateWithRBITable[Climate %in% c("ATA", "GSTA", "NONGSTA"), ClimateName:="Temperature"]
climateWithRBITable[Climate %in% c("ACMIA", "GSCMIA", "NONGSCMIA"), ClimateName:="CMI"]
climateWithRBITable[Climate %in% c("ACO2A", "GSCO2A", "NONGSCO2A"), ClimateName:="CO2"]
climateWithRBITable[, ClimateName:=factor(ClimateName, levels = c("Temperature", "CMI", "CO2"))]

rm(interactionTable, maineffectTable, speciesdata)
maineffectTable <- unique(climateWithRBITable[maineffect != 0,.(Species, Climate, Season, ClimateName,
                                                 climateEffect = maineffect, maineffect_SE, RBI = meanRBI)],
                          by = c("Species", "Climate"))


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


