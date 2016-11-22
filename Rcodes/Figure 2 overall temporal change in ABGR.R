rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "Results", "finalYearModels.RData"))
rm(speciesData, indispecies)
##### for overall temporal trends and its dependency on DBH and RBI
studySpecies <- c("All", "JP", "TA", "BS", "Other")
for(indispecies in studySpecies){
  speciesData <- finalAnalysesData[[indispecies]]
  overallTrendOutput <- data.table(Year = seq(min(speciesData$Year), 
                                              max(speciesData$Year),
                                              length = 100))
  overallTrendOutput[,':='(Species = indispecies,
                           Yearctd = Year-mean(speciesData$Year),
                           RBIctd = 0,
                           logDBHctd = 0,
                           logIntraHctd = 0,
                           logInterHctd = 0,
                           logSBctd = 0,
                           logSActd = 0)]
  themodel <- finalModels[[indispecies]]
  # for overall trend
  fittedvalues <- predict(themodel, newdata = overallTrendOutput, level = 0, se.fit = TRUE)
  overallTrendOutput$PredictedABGR <- exp(fittedvalues$fit)
  overallTrendOutput$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  overallTrendOutput$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  rm(fittedvalues)
  overallTrendOutput <- overallTrendOutput[,.(Species, Direction = "Overall", Year, 
                                                PredictedABGR,
                                                PredictedABGR_Lower, PredictedABGR_Upper)]
  if(indispecies == "All"){
    allFigureData <- overallTrendOutput
  } else {
    allFigureData <- rbind(allFigureData, overallTrendOutput)
  }
  rm(overallTrendOutput)
}
rm(indispecies)
allFigureData[,Species:=factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                                     labels = c("All species", "Jack pine",
                                                "Trembling aspen", "Black spruce",
                                                "Other species"))]
allFigureData <- allFigureData[Species != "Trembling aspen",]
allFigureData[Species %in% c("All species", "Jack pine"), 
              yfacet:=factor("a", levels = letters[1:2])]
allFigureData[Species %in% c("Black spruce", "Other species"), 
              yfacet:=factor("b", levels = letters[1:2])]
allFigureData[Species %in% c("All species", "Black spruce"), 
              xfacet := "a"]
allFigureData[Species %in% c("Jack pine", "Other species"), 
              xfacet := "b"]


texts <- data.table(xfacet = c(rep("a", 2), rep("b", 2)),
                    yfacet = rep(c("a", "b"), 2),
                    x = 1986, y = rep(c(1.1, 0.6), 2),
                    Species = c("a", "c",
                                "b", "d"))


Figure2 <- ggplot(data = allFigureData[!(Species %in% c(" ", "Trembling aspen"))], aes(x = Year, y = PredictedABGR))+
  facet_grid(yfacet~xfacet, scales = "free_y")+
  # for a overall trends
  geom_ribbon(aes(x = Year, ymin = PredictedABGR_Lower, ymax = PredictedABGR_Upper),
              fill = "blue", alpha = 0.1)+
  geom_line(aes(x = Year, y = PredictedABGR), size = 1, col = "blue")+
  geom_text(data = texts,
             aes(x = x, y = y, label = Species), hjust = 0, size = 10)+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  scale_y_continuous(name = expression(paste("Aboveground biomass growth rate (Kg ", year^{-1}, ")")))+
  geom_segment(aes(x=-Inf, xend=-Inf, y=Inf, yend=-Inf), size = 1.5)+
  geom_segment(aes(y=-Inf, yend=-Inf, x=Inf, xend=-Inf), size = 1.5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        # axis.line.x = element_line(size = 1, colour = "black"),
        # axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text = element_blank())


ggsave(file = file.path(workPath, "TablesFigures", "Figure 2_temporal trends in ABGR.png"), Figure2,
       width = 10, height = 9.2)
