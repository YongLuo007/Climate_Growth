rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "Results", "bestYearModel.RData"))
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table

ThreeDGrowthvsYearandDom <- data.table(Species = character(), RBI = numeric(),
                                       Year = numeric(), PredictedBGR = numeric(),
                                       PredictedBGR_Upper = numeric(), PredictedBGR_Lower = numeric(),
                                       Main = numeric())
for(indispecies in c("All", "JP", "BS", "TA", "Other")){
  themodel <- allbestmodels[[paste(indispecies)]]
  if(indispecies == "All"){
    newdata <- data.table(expand.grid(Year = seq(min(analysesData$Year), 
                                                 max(analysesData$Year),
                                                 length = 100),
                                      RBI = c(seq(0, 100, by = 10),
                                              mean(analysesData$RBI))))
    newdata[,':='(Yearctd = Year-mean(analysesData$Year),
                RBIctd = RBI - mean(analysesData$RBI),
                logDBHctd = log(mean(analysesData$IniDBH))-
                  mean(log(analysesData$IniDBH)),
                logHctd = log(mean(analysesData$Hegyi))-
                  mean(log(analysesData$Hegyi)),
                Species = indispecies)]
    newdata2 <- data.table::copy(analysesData)
  } else if (indispecies == "Other"){
    newdata <- data.table(expand.grid(Year = seq(min(analysesData[!(Species %in% c("JP", "TA", "BS")),]$Year), 
                                                 max(analysesData[!(Species %in% c("JP", "TA", "BS")),]$Year),
                                                 length = 100),
                                      RBI = c(seq(0, 100, by = 10),
                                              mean(analysesData[!(Species %in% c("JP", "TA", "BS")),]$RBI))))
    newdata[,':='(Yearctd = Year-mean(analysesData[!(Species %in% c("JP", "TA", "BS")),]$Year),
                  RBIctd = RBI - mean(analysesData[!(Species %in% c("JP", "TA", "BS")),]$RBI),
                  logDBHctd = log(mean(analysesData[!(Species %in% c("JP", "TA", "BS")),]$IniDBH))-
                    mean(log(analysesData[!(Species %in% c("JP", "TA", "BS")),]$IniDBH)),
                  logHctd = log(mean(analysesData[!(Species %in% c("JP", "TA", "BS")),]$Hegyi))-
                    mean(log(analysesData[!(Species %in% c("JP", "TA", "BS")),]$Hegyi)),
                  Species = indispecies)]
    newdata2 <- analysesData[!(Species %in% c("JP", "TA", "BS")), ]
    
  } else {
    newdata <- data.table(expand.grid(Year = seq(min(analysesData[Species == indispecies,]$Year), 
                                                 max(analysesData[Species == indispecies,]$Year),
                                                 length = 100),
                                      RBI = c(seq(0, 100, by = 10),
                                              mean(analysesData[Species == indispecies,]$RBI))))
    newdata[,':='(Yearctd = Year-mean(analysesData[Species == indispecies,]$Year),
                  RBIctd = RBI - mean(analysesData[Species == indispecies,]$RBI),
                  logDBHctd = log(mean(analysesData[Species == indispecies,]$IniDBH))-
                    mean(log(analysesData[Species == indispecies,]$IniDBH)),
                  logHctd = log(mean(analysesData[Species == indispecies,]$Hegyi))-
                    mean(log(analysesData[Species == indispecies,]$Hegyi)),
                  Species = indispecies)]
    newdata2 <- analysesData[Species == indispecies, ]
  }
  newdata2[,':='(Yearctd = Year-mean(Year), RBIctd = 0, 
                 logDBHctd = log(mean(IniDBH))-mean(log(IniDBH)),
                 logHctd = log(mean(Hegyi))-mean(log(Hegyi)))]
  
  reducedFomu <- as.formula(paste(as.character(formula(summary(themodel)$term))[c(2, 1, 3)], collapse = " "))
  fittedvalues <- predict(themodel, newdata = newdata, level = 0, se.fit = TRUE)
  # data.frame(summary(themodel)$tTable)$Value[1]
  newdata2[, predictedBGR := predict(themodel, newdata = newdata2, 
                                     level = 0)]
  addRandom <- FALSE
  if(addRandom){
    plotrandom <- ranef(themodel)[[1]]
    names(plotrandom)[1] <- "PlotIntercept"
    plotrandom$PlotID <- row.names(plotrandom)
    plotrandom <- data.table(plotrandom)[,.(PlotID, PlotIntercept)]
    treerandom <- ranef(themodel)[[2]]
    names(treerandom)[1] <- "TreeIntercept"
    treerandom$uniTreeID <- row.names(treerandom)
    treerandom <- data.table(treerandom)[,.(uniTreeID, TreeIntercept)]
    treerandom[, uniTreeID:=unlist(lapply(strsplit(uniTreeID, "/", fixed = TRUE), function(x){x[2]}))]
    newdata2 <- setkey(newdata2, PlotID)[setkey(plotrandom, PlotID),
                                         nomatch = 0]
    newdata2 <- setkey(newdata2, uniTreeID)[setkey(treerandom, uniTreeID),
                                            nomatch = 0]
    newdata2[,predictedBGR:=exp(predictedBGR+PlotIntercept+TreeIntercept)]
    set(newdata2, ,c("PlotIntercept", "TreeIntercept"), NULL)
  } else {
    newdata2[,predictedBGR:=exp(predictedBGR)]
  }
  newdata2 <- newdata2[, .(totalPredictedBGR = sum(predictedBGR), meanPredBGR = mean(predictedBGR),
                           Species = indispecies), by = c("PlotID", "Year")]
  newdata$PredictedBGR <- exp(fittedvalues$fit)
  newdata$PredictedBGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  newdata$PredictedBGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  newdata[,Main:=0]
  if(indispecies == "All"){
    newdata[RBI == mean(analysesData$RBI), Main := 1]
  } else if(indispecies == "Other"){
    newdata[RBI == mean(analysesData[!(Species %in% c("JP", "TA", "BS")),]$RBI), Main := 1]
  } else {
    newdata[RBI == mean(analysesData[Species == indispecies,]$RBI), Main := 1]
  }
  
  ThreeDGrowthvsYearandDom <- rbind(ThreeDGrowthvsYearandDom,
                                    newdata[,.(Species, Year, RBI, PredictedBGR, 
                                               PredictedBGR_Upper, PredictedBGR_Lower, Main)])
  rm(themodel, newdata, fittedvalues, reducedFomu)
  if(indispecies == "All"){
    PlotTrendData <- newdata2
  } else {
    PlotTrendData <- rbind(PlotTrendData, newdata2)
  }
}
PlotTrendData[,Species:=factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                               labels = c("All species", "Jack pine",
                                          "Trembling aspen", "Black spruce",
                                          "Other species"))]
a <- ggplot(data=PlotTrendData, aes(x = Year, y = totalPredictedBGR))+
  geom_line(aes(group = PlotID), col = "gray")+
  facet_grid(Species~., scales = "free_y")

b <- ggplot(data=PlotTrendData, aes(x = Year, y = meanPredBGR))+
  geom_line(aes(group = PlotID), col = "gray")+
  facet_grid(Species~.)


ThreeDGrowthvsYearandDom$Species <- factor(ThreeDGrowthvsYearandDom$Species,
                                           levels = c("All", "JP", "TA", "BS"),
                                           labels = c("All species","Jack pine", "Trembling aspen", "Black spruce"))
OtherVariable$Species <- factor(OtherVariable$Species,
                                levels = c("All", "JP", "TA", "BS"),
                                labels = c("All species","Jack pine", "Trembling aspen", "Black spruce"))
OtherVariable[, ':='(y1 = 5, y2 = 4, Year = 1993, DBH1 = paste(DBH, "cm"))]
MainTrend <- ThreeDGrowthvsYearandDom[Main == 1,][,linetype := 1]
MainTrend[Species %in% c("Black spruce"), linetype := 2]
MainTrend[,linetype:=factor(linetype, levels=c(1, 2))]
Fig2_left <- ggplot(data = ThreeDGrowthvsYearandDom[Main == 0], aes(x = Year, y = PredictedBGR))+
  geom_line(aes(group = RBI, col = RBI))+
  scale_colour_continuous(name = "RBI", low = "#FF0000", high = "#00FF00", breaks = seq(0, 100, by = 20))+
  geom_text(data = data.frame(Year = rep(1985, 4), y = c(1.5, 1.5, 1.5, 1.5), label = c("a", "b", "c", "d"),
                              Species = c("All species", "Jack pine", "Trembling aspen", "Black spruce")),
            aes(x = Year, y = y, label = label), size = 8)+
  geom_line(data = MainTrend, aes(x = Year, y = PredictedBGR, linetype = linetype),
            colour = "black", size = 1)+
  facet_grid(Species~.)+
  scale_linetype(guide = "none")+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  scale_y_continuous(name = expression(atop("Aboveground biomass growth rate", paste("(Kg ", year^{-1}, ")"))),
                     limits = c(0, 1.5), breaks = seq(0, 1.5, 0.3))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1)+
  # annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.background = element_blank(),
        # strip.text = element_blank(),
        legend.position = c(0.1, 0.8),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

for(indispecies in c("JP", "TA", "BS")){
  a <- density(analysesData[Species == indispecies, ]$Dominance_indiBiomass,
               adjust = 1.5, from = 0, to = 100)
  densitydata <- data.table(Species = indispecies, Dominance = a$x, Density = a$y)
  if(indispecies == "JP"){
    alldensitydata <- densitydata
  } else{
    alldensitydata <- rbind(alldensitydata, densitydata)
  }
}
rm(a, indispecies, densitydata)
alldensitydata[, Species := factor(Species, levels=c("JP", "TA", "BS"),
                                   labels = c("Jack pine", "Trembling aspen", "Black spruce"))]
meanDomData <- analysesData[,.(minDom = mean(Dominance_indiBiomass)), by = Species][
  ,':='(Species = factor(Species, levels=c("JP", "TA", "BS"),
                         labels = c("Jack pine", "Trembling aspen", "Black spruce")),
        Density0 = -0.03, Density10 = 0.03)]
Fig2_right <- ggplot(data = alldensitydata[,':='(Density0 = -Density/2, Density10 = Density/2)],
                     aes(x = Dominance, y = Density0))+
  geom_segment(aes(xend = Dominance, yend = Density10, col = Dominance))+
  geom_segment(data = meanDomData, aes(x = minDom, xend = minDom, y = Density0, yend = Density10),
               colour = "gray", size = 1)+
  facet_grid(Species~., scales = "free_y")+
  geom_text(data = data.frame(Year = rep(100, 3), y = c(-0.035, -0.035, -0.035), label = c("b", "d", "f"),
                              Species = c("Jack pine", "Trembling aspen", "Black spruce")),
            aes(x = Year, y = y, label = label), size = 8)+
  coord_flip()+
  scale_colour_continuous(low = "#FF0000", high = "#00FF00", breaks = seq(0, 100, by = 20))+
  scale_y_continuous(name = "Density", limits = c(-0.04, 0.04),
                     breaks = seq(-0.04, 0.04, by = 0.04))+
  scale_x_continuous(name = "RBI")+
  theme_bw()+
  # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        # axis.ticks.y = element_blank(),
        strip.background = element_rect(colour = "white", fill = "gray"),
        strip.text = element_text(size = 12),
        legend.position = "none")

Fig2_left_Grob <- ggplotGrob(Fig2_left)
Fig2_right_Grob <- ggplotGrob(Fig2_right)

Fig2_right_Grob$heights <- Fig2_left_Grob$heights

dev(4)
clearPlot()
plotlayout <- rbind(c(1, 1, 1, 2))
c <- grid.arrange(Fig2_left_Grob, Fig2_right_Grob,
                  layout_matrix = plotlayout)

ggsave(file = file.path(workPath, "TablesFigures", "Fig2.png"), c,
       width = 10, height = 10)



