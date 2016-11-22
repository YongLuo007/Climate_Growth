rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "Results", "finalYearModels.RData"))
rm(speciesData, indispecies)
##### for overall temporal trends and its dependency on DBH and RBI
for(indispecies in studySpecies){
  speciesData <- finalAnalysesData[[indispecies]]
  changewithIntraH <- data.table(expand.grid(Species = indispecies,
                                             Direction = "changewithIntraH",
                                             Year = seq(min(speciesData$Year), 
                                                        max(speciesData$Year),
                                                        length = 100), 
                                             IntraH = c(min(speciesData$IntraH),
                                                        exp(mean(log(speciesData$IntraH+1)))-1,
                                                        max(speciesData$IntraH)),
                                             stringsAsFactors = FALSE))
  changewithIntraH[,':='(Yearctd = Year-mean(speciesData$Year),
                         logDBHctd = 0,
                         logIntraHctd = log(IntraH+1)-mean(log(speciesData$IntraH+1)),
                         logInterHctd = 0,
                         logSActd = 0)]
  changewithIntraH[IntraH == min(speciesData$IntraH), LineCol:="Minimum"]
  changewithIntraH[IntraH == exp(mean(log(speciesData$IntraH+1)))-1, LineCol:="Mean"]
  changewithIntraH[IntraH == max(speciesData$IntraH), LineCol:="Maximum"]
  
  changewithInterH <- data.table(expand.grid(Species = indispecies,
                                             Direction = "changewithInterH",
                                             Year = seq(min(speciesData$Year), 
                                                        max(speciesData$Year),
                                                        length = 100), 
                                             InterH = c(min(speciesData$InterH),
                                                        exp(mean(log(speciesData$InterH+1)))-1,
                                                        max(speciesData$InterH)),
                                             stringsAsFactors = FALSE))
  changewithInterH[,':='(Yearctd = Year-mean(speciesData$Year),
                         logDBHctd = 0,
                         logIntraHctd = 0,
                         logInterHctd = log(InterH+1)-mean(log(speciesData$InterH+1)),
                         logSActd = 0)]
  changewithInterH[InterH == min(speciesData$InterH), LineCol:="Minimum"]
  changewithInterH[InterH == exp(mean(log(speciesData$InterH+1)))-1, LineCol:="Mean"]
  changewithInterH[InterH == max(speciesData$InterH), LineCol:="Maximum"]
  
  
  changewithDBH <- data.table(expand.grid(Species = indispecies,
                                          Direction = "changewithDBH",
                                          Year = seq(min(speciesData$Year), 
                                                     max(speciesData$Year),
                                                     length = 100), 
                                          DBH = c(min(speciesData$IniDBH),
                                                  exp(mean(log(speciesData$IniDBH))),
                                                  max(speciesData$IniDBH)),
                                          stringsAsFactors = FALSE))
  changewithDBH[,':='(Yearctd = Year-mean(speciesData$Year),
                      logDBHctd = log(DBH)-mean(log(speciesData$IniDBH)),
                      logIntraHctd = 0,
                      logInterHctd = 0,
                      logSActd = 0)]
  changewithDBH[DBH == min(speciesData$IniDBH), LineCol:="Minimum"]
  changewithDBH[DBH == exp(mean(log(speciesData$IniDBH))), LineCol:="Mean"]
  changewithDBH[DBH == max(speciesData$IniDBH), LineCol:="Maximum"]
  
  
  changewithSA <- data.table(expand.grid(Species = indispecies,
                                         Direction = "changewithSA",
                                         Year = seq(min(speciesData$Year), 
                                                    max(speciesData$Year),
                                                    length = 100), 
                                         SA = c(min(speciesData$SA),
                                                exp(mean(log(speciesData$SA))),
                                                max(speciesData$SA)),
                                         stringsAsFactors = FALSE))
  changewithSA[,':='(Yearctd = Year-mean(speciesData$Year),
                     logDBHctd = 0,
                     logIntraHctd = 0,
                     logInterHctd = 0,
                     logSActd = log(SA)-mean(log(speciesData$SA)))]
  changewithSA[SA == min(speciesData$SA), LineCol:="Minimum"]
  changewithSA[SA == exp(mean(log(speciesData$SA))), LineCol:="Mean"]
  changewithSA[SA == max(speciesData$SA), LineCol:="Maximum"]
  
  themodel <- finalModels[[indispecies]]
  fittedvalues <- predict(themodel, newdata = changewithIntraH, level = 0, se.fit = TRUE)
  changewithIntraH$PredictedABGR <- exp(fittedvalues$fit)
  changewithIntraH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  changewithIntraH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  rm(fittedvalues)
  
  fittedvalues <- predict(themodel, newdata = changewithInterH, level = 0, se.fit = TRUE)
  changewithInterH$PredictedABGR <- exp(fittedvalues$fit)
  changewithInterH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  changewithInterH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  rm(fittedvalues)
  
  fittedvalues <- predict(themodel, newdata = changewithDBH, level = 0, se.fit = TRUE)
  changewithDBH$PredictedABGR <- exp(fittedvalues$fit)
  changewithDBH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  changewithDBH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  rm(fittedvalues)
  
  fittedvalues <- predict(themodel, newdata = changewithSA, level = 0, se.fit = TRUE)
  changewithSA$PredictedABGR <- exp(fittedvalues$fit)
  changewithSA$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  changewithSA$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  rm(fittedvalues)
  
  speciesOutput <- rbind(changewithIntraH[,.(Species, Direction, Year, LineCol, PredictedABGR, PredictedABGR_Lower, PredictedABGR_Upper)],
                         changewithInterH[,.(Species, Direction, Year, LineCol, PredictedABGR, PredictedABGR_Lower, PredictedABGR_Upper)],
                         changewithDBH[,.(Species, Direction, Year, LineCol, PredictedABGR, PredictedABGR_Lower, PredictedABGR_Upper)],
                         changewithSA[,.(Species, Direction, Year, LineCol, PredictedABGR, PredictedABGR_Lower, PredictedABGR_Upper)])
  
  if(indispecies == "All species"){
    allFigureData <- speciesOutput
  } else {
    allFigureData <- rbind(allFigureData, speciesOutput)
  }
}

allFigureData[,':='(Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen",
                                                "Black spruce", "Other species")),
                    LineCol = factor(LineCol, levels = c("Maximum", "Mean", "Minimum")),
                    Direction = factor(Direction, levels = c("changewithIntraH", "changewithInterH",
                                                             "changewithDBH", "changewithSA")))]

allFigureData1 <- allFigureData[Species != "All species" |
                                  !(Direction %in% c("changewithSA")),]
allFigureData1 <- allFigureData1[Species != "Jack pine" | 
                                   !(Direction %in% c("changewithInterH", "changewithDBH")),]

allFigureData1 <- allFigureData1[Species != "Trembling aspen" | 
                                   !(Direction %in% c("changewithIntraH", "changewithInterH")),]

allFigureData1 <- allFigureData1[Species != "Black spruce" | 
                                   !(Direction %in% c("changewithSA")),]
allFigureData1 <- allFigureData1[Species != "Other species" |
                                   !(Direction %in% c("changewithIntraH", "changewithInterH",
                                                      "changewithDBH")),]
segmenttable <- data.table(Species = "All species", Direction = factor(c("changewithIntraH", "changewithInterH",
                                                                         "changewithDBH", "changewithSA")),
                           x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
texttable <- data.table(Species = "All species", Direction = factor(c("changewithIntraH", "changewithInterH",
                                                                      "changewithDBH", "changewithSA")),
                        x = 1989, y = Inf, texts = letters[1:4])
texttable <- rbind(texttable, data.table(Species = c("Jack pine", "Trembling aspen",
                                                     "Black spruce", "Other species"),
                                         Direction = "changewithIntraH", x = 1989, y = Inf, texts = " "))
texttable[,':='(Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen",
                                                     "Black spruce", "Other species")),
                Direction = factor(Direction, levels = c("changewithIntraH", "changewithInterH",
                                                         "changewithDBH", "changewithSA")))]
segmenttable2 <- data.table(expand.grid(Species = factor(c("All species", "Jack pine", "Trembling aspen",
                                                           "Black spruce", "Other species")),
                                        Direction = factor(c("changewithIntraH", "changewithInterH",
                                                             "changewithDBH", "changewithSA"))))
segmenttable <- rbind(segmenttable, segmenttable2[,':='(x = -Inf, xend = Inf, y = -Inf, yend = -Inf)])
segmenttable[,':='(Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen",
                                                        "Black spruce", "Other species")),
                   Direction = factor(Direction, levels = c("changewithIntraH", "changewithInterH",
                                                            "changewithDBH", "changewithSA")))]

figure <- ggplot(data = allFigureData1, aes(x = Year, y = PredictedABGR))+
  facet_grid(Direction~Species, scales = "free_y")+
  geom_line(aes(group = LineCol, col = LineCol), size = 1)+
  # geom_line(aes(group = LineCol, y = PredictedABGR_Lower, col = LineCol), linetype = 2)+
  # geom_line(aes(group = LineCol, y = PredictedABGR_Upper, col = LineCol), linetype = 2)+
  scale_y_continuous(name = expression(paste("Aboveground biomass growth rate (Kg ", year^{-1}, ")")))+
  scale_x_continuous(name = "Year", breaks = seq(1990, 2010, by = 5))+
  scale_color_manual(name = "Competition intensity", values = c("blue", "black", "magenta"),
                     labels = c("Strong", "Medium strong", "Weak"))+
  guides(col = guide_legend(title.position = "top",
                            direction = "horizontal"))+
  geom_segment(data = segmenttable, aes(x = x, xend = xend, y = y, yend = yend), size = 1)+
  geom_text(data = texttable, aes(x = x, y = y, label = texts), vjust = 1.5, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # axis.line.x = element_line(size = 1, colour = "black"),
        # axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour = "white"),
        strip.text.x = element_text(size = 15),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.80, 0.85),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))


dev(4)
clearPlot()

ggsave(file = file.path(workPath, "TablesFigures", "Figure 3_temporal trends by competition.png"),
       figure,  width = 12, height = 8)






