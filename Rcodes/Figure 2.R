rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "Results", "bestYearModel.RData"))
rm(speciesData, indispecies)
##### for overall temporal trends and its dependency on DBH and RBI

for(indispecies in studySpecies){
  if(indispecies == "All"){
    speciesData <- data.table::copy(analysesData)
  } else if(indispecies == "Other"){
    speciesData <- analysesData[!(Species %in% c("JP", "TA", "BS")),]
  } else {
    speciesData <- analysesData[Species == indispecies,]
  }
  overallTrendOutput <- data.table(Year = seq(min(speciesData$Year), 
                                              max(speciesData$Year),
                                              length = 100))
  TrendWithRBI <- data.table(expand.grid(Year = seq(min(speciesData$Year), 
                                                    max(speciesData$Year),
                                                    length = 100),
                                         RBI = seq(0, 100, length = 10)))
  TrendWithDBH <- data.table(expand.grid(Year = seq(min(speciesData$Year), 
                                                    max(speciesData$Year),
                                                    length = 100),
                                         DBH = seq(min(speciesData$IniDBH),
                                                   max(speciesData$IniDBH),
                                                   length = 10)))
  overallTrendOutput[,':='(Species = indispecies,
                           Yearctd = Year-mean(speciesData$Year),
                           RBIctd = 0,
                           logDBHctd = log(mean(speciesData$IniDBH))-
                             mean(log(speciesData$IniDBH)),
                           logHctd = log(mean(speciesData$Hegyi))-
                             mean(log(speciesData$Hegyi)))]
  TrendWithRBI[,':='(Species = indispecies,
                     Yearctd = Year-mean(speciesData$Year),
                     RBIctd = RBI-mean(speciesData$RBI),
                     logDBHctd = log(mean(speciesData$IniDBH))-
                       mean(log(speciesData$IniDBH)),
                     logHctd = log(mean(speciesData$Hegyi))-
                       mean(log(speciesData$Hegyi)))]
  TrendWithDBH[,':='(Species = indispecies,
                     Yearctd = Year-mean(speciesData$Year),
                     RBIctd = 0,
                     logDBHctd = log(DBH)-
                       mean(log(speciesData$IniDBH)),
                     logHctd = log(mean(speciesData$Hegyi))-
                       mean(log(speciesData$Hegyi)))]
  themodel <- allbestmodels[[indispecies]]
  # for overall trend
  fittedvalues <- predict(themodel, newdata = overallTrendOutput, level = 0, se.fit = TRUE)
  overallTrendOutput$PredictedABGR <- exp(fittedvalues$fit)
  overallTrendOutput$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  overallTrendOutput$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  rm(fittedvalues)
  # for trend dependent on RBI
  TrendWithRBI[,PredictedABGR := exp(predict(themodel, newdata = TrendWithRBI,
                                             level = 0, se.fit = FALSE))]
  # for trend dependent on DBH
  TrendWithDBH[,PredictedABGR := exp(predict(themodel, newdata = TrendWithDBH,
                                             level = 0, se.fit = FALSE))]
  onespeciesData <- rbind(overallTrendOutput[,.(Species, Direction = "Overall", Year, RBI = 0,
                                                DBH = 0, PredictedABGR,
                                                PredictedABGR_Lower, PredictedABGR_Upper)],
                          TrendWithRBI[,.(Species, Direction = "withRBI", Year, RBI,
                                          DBH = 0, PredictedABGR, PredictedABGR_Lower = 0,
                                          PredictedABGR_Upper = 0)],
                          TrendWithDBH[,.(Species, Direction = "withDBH", Year, RBI = 0,
                                          DBH, PredictedABGR, PredictedABGR_Lower = 0,
                                          PredictedABGR_Upper = 0)])
  if(indispecies == "All"){
    allFigureData <- onespeciesData
  } else {
    allFigureData <- rbind(allFigureData, onespeciesData)
  }
  rm(overallTrendOutput, TrendWithDBH, TrendWithRBI)
}
rm(indispecies)

allFigureData[,Species:=factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                                     labels = c("All species", "Jack pine",
                                                "Trembling aspen", "Black spruce",
                                                "Other species"))]
allFigureData[, Direction:=factor(Direction, levels = c("Overall", "withRBI", "withDBH"),
                                  labels = c("Overall trend", "Change with RBI", "Change with DBH"))]

Figure2 <- ggplot(data = allFigureData[Direction %in% c("Overall trend", "Change with RBI"),], aes(x = Year, y = PredictedABGR))+
  # for a overall trends
  geom_ribbon(data = allFigureData[Direction == "Overall trend",],
              aes(x = Year, ymin = PredictedABGR_Lower, ymax = PredictedABGR_Upper),
              fill = "blue", alpha = 0.1)+
  geom_line(data = allFigureData[Direction == "Overall trend",],
            aes(x = Year, y = PredictedABGR), size = 1, col = "blue")+
  # for b change with RBI
  geom_line(data = allFigureData[Direction == "Change with RBI" & Species != "Other species",],
            aes(x = Year, y = PredictedABGR, group = RBI, col = RBI),
            size = 1)+
  # accessories
  geom_segment(aes(x=-Inf, xend=-Inf, y=Inf, yend=-Inf), size = 1.5)+
  geom_text(data = data.frame(Year = rep(1986, 2), y = 1.45, label = c("a", "b"), 
                              Species = "All species", Direction = c("Overall trend", "Change with RBI")),
            aes(x = Year, y = y, label = label), size = 10)+
  
  scale_colour_continuous(name = "RBI", low = "#FF0000", high = "#00FF00", breaks = c(0, 50, 100))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  scale_y_continuous(name = expression(paste("Aboveground biomass growth rate (Kg  ", year^{-1}, ")")),
                     limits = c(0, 1.52), breaks = seq(0, 1.5, 0.5))+
  
  facet_grid(Species~Direction)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.background = element_rect(colour = "white"),
        strip.text.y = element_text(size = 16),
        # strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.75, 0.1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))


ggsave(file = file.path(workPath, "TablesFigures", "Figure 2_temporal trends in ABGR.png"), Figure2,
       width = 10, height = 10)
