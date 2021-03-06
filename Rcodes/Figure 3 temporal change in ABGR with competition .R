rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "finalYearModels.RData"))
##### for overall temporal trends and its dependency on DBH and RBI
for(i in 1:length(allFixedCoeff)){
  indicoeff <- allFixedCoeff[[i]]
  indicoeff[,Species:=names(allFixedCoeff[i])]
  if(i == 1){
    allcoeff <- indicoeff
  } else {
    allcoeff <- rbind(allcoeff, indicoeff)
  }
}
allcoeff[, variable:= unlist(lapply(lapply(lapply((rn), function(x) unlist(strsplit(x, ":"))), function(y) sort(y)),
                                    function(z) paste(z, collapse = ":")))]
# unique(allcoeff$variable)
allcoeff[variable == "logDBHctd:Yearctd", Direction:="changewithDBH"]
allcoeff[variable == "logSActd:Yearctc", Direction:="changewithSA"]
allcoeff[variable == "logIntraHctd:Yearctd", Direction:="changewithIntraH"]
allcoeff[variable == "logInterHctd:Yearctd", Direction:="changewithInterH"]
output <- data.table(Species = character(), Direction = character(), Year = numeric(),
                     CompetitionIntensity = character(), PredictedABGR = numeric(),
                     PredictedABGR_Lower = numeric(), PredictedABGR_Upper = numeric())
for(indispecies in studySpecies){
  speciesData <- analysesData[DataType == indispecies,]
  speciecoeff <- allcoeff[Species == indispecies, ]
  bestFormu <- bestFomula[[indispecies]]
  themodel <- bestModels[[indispecies]]
  if(nrow(speciecoeff[rn == "Yearctd",]) == 1){
    mainTrends <- data.table(Species = indispecies, 
                             Direction = "mainTrend",
                             Year = seq(min(speciesData$Year), 
                                        max(speciesData$Year),
                                        length = 100))
    mainTrends[,':='(Yearctd = Year - mean(speciesData$Year),
                     logDBHctd = 0,
                     logSActd = 0,
                     logIntraHctd = 0,
                     logInterHctd = 0)]
    fittedvalues <- predict(themodel, newdata = mainTrends, level = 0, se.fit = TRUE)
    mainTrends$PredictedABGR <- exp(fittedvalues$fit)
    mainTrends$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    mainTrends$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    output <- rbind(output, mainTrends[,.(Species, Direction, Year, 
                                          CompetitionIntensity = "Medium",
                                          PredictedABGR, PredictedABGR_Lower,
                                          PredictedABGR_Upper)])
    rm(fittedvalues)
  }
  if(nrow(speciecoeff[Direction == "changewithIntraH",]) == 1){
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
    changewithIntraH[IntraH == min(speciesData$IntraH), CompetitionIntensity:="Weak"]
    changewithIntraH[IntraH == exp(mean(log(speciesData$IntraH+1)))-1, CompetitionIntensity:="Medium"]
    changewithIntraH[IntraH == max(speciesData$IntraH), CompetitionIntensity:="Strong"]
    fittedvalues <- predict(themodel, newdata = changewithIntraH, level = 0, se.fit = TRUE)
    changewithIntraH$PredictedABGR <- exp(fittedvalues$fit)
    changewithIntraH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithIntraH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    rm(fittedvalues)
    output <- rbind(output, changewithIntraH[,.(Species, Direction, Year, 
                                                CompetitionIntensity,
                                                PredictedABGR, PredictedABGR_Lower, PredictedABGR_Upper)])
  }
  if(nrow(speciecoeff[Direction == "changewithInterH",]) == 1){
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
    changewithInterH[InterH == min(speciesData$InterH), CompetitionIntensity:="Weak"]
    changewithInterH[InterH == exp(mean(log(speciesData$InterH+1)))-1, CompetitionIntensity:="Medium"]
    changewithInterH[InterH == max(speciesData$InterH), CompetitionIntensity:="Strong"]
    fittedvalues <- predict(themodel, newdata = changewithInterH, level = 0, se.fit = TRUE)
    changewithInterH$PredictedABGR <- exp(fittedvalues$fit)
    changewithInterH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithInterH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    rm(fittedvalues)
    output <- rbind(output, changewithInterH[,.(Species, Direction, Year, 
                                                CompetitionIntensity,
                                                PredictedABGR, PredictedABGR_Lower,
                                                PredictedABGR_Upper)])
  }
  if(nrow(speciecoeff[Direction == "changewithDBH",]) == 1){
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
    changewithDBH[DBH == min(speciesData$IniDBH), CompetitionIntensity:="Weak"]
    changewithDBH[DBH == exp(mean(log(speciesData$IniDBH))), CompetitionIntensity:="Medium"]
    changewithDBH[DBH == max(speciesData$IniDBH), CompetitionIntensity:="Strong"]
    fittedvalues <- predict(themodel, newdata = changewithDBH, level = 0, se.fit = TRUE)
    changewithDBH$PredictedABGR <- exp(fittedvalues$fit)
    changewithDBH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithDBH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    rm(fittedvalues)
    output <- rbind(output, changewithDBH[,.(Species, Direction, Year, 
                                                CompetitionIntensity,
                                                PredictedABGR, PredictedABGR_Lower,
                                                PredictedABGR_Upper)])
    
  }
  if(nrow(speciecoeff[Direction == "changewithSA",]) == 1){
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
    changewithSA[SA == min(speciesData$SA), CompetitionIntensity:="Weak"]
    changewithSA[SA == exp(mean(log(speciesData$SA))), CompetitionIntensity:="Medium"]
    changewithSA[SA == max(speciesData$SA), CompetitionIntensity:="Strong"]
    fittedvalues <- predict(themodel, newdata = changewithSA, level = 0, se.fit = TRUE)
    changewithSA$PredictedABGR <- exp(fittedvalues$fit)
    changewithSA$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithSA$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    rm(fittedvalues)
    output <- rbind(output, changewithSA[,.(Species, Direction, Year, 
                                             CompetitionIntensity,
                                             PredictedABGR, PredictedABGR_Lower,
                                             PredictedABGR_Upper)])
    
  }
}
output <- output[Species != "All species",]
output[Species == "Other species", Species:="Minor species"]
fourspecies <- c("Jack pine", "Trembling aspen", "Black spruce", "Minor species")
output[,':='(Species = factor(Species, fourspecies))]
plotFigure <- "CompetitionOnly"
if(plotFigure == "CompetitionOnly"){
  directions <- c("mainTrend", "changewithIntraH", "changewithInterH")
} else if (plotFigure == "OntogenyOnly") {
  directions <- c("mainTrend", "changewithDBH", "changewithSA")
} else {
  directions <- c("mainTrend", "changewithIntraH", "changewithInterH", "changewithDBH", "changewithSA")
}
allFigureData <- data.table::copy(output)[Direction %in% directions, ]
allFigureData[, ':='(Direction = factor(Direction, levels = directions),
                     CompetitionIntensity = factor(CompetitionIntensity, 
                                                   levels = c("Strong", "Medium", "Weak")))]

segmenttable <- data.table(Species = "Jack pine", 
                           Direction = directions,
                           x = -Inf, xend = -Inf,yend = Inf)

segmenttable2 <- data.table(expand.grid(Species = factor(fourspecies),
                                        Direction = factor(directions, levels = directions)))
segmenttable2[,':='(x = -Inf, xend = Inf)]
segmenttable2[Direction == "mainTrend", ':='(yend = 0.134)]
segmenttable2[Direction == "changewithIntraH", ':='(yend = 0.001)]
segmenttable2[Direction == "changewithInterH", ':='(yend = 0.026)]

segmenttable <- rbind(segmenttable, segmenttable2)
segmenttable[Direction == "mainTrend", ':='(y = 0.134)]
segmenttable[Direction == "changewithIntraH", ':='(y = 0.001)]
segmenttable[Direction == "changewithInterH", ':='(y = 0.026)]


segmenttable[,':='(Species = factor(Species, 
                                    levels = fourspecies),
                   Direction = factor(Direction, 
                                      levels = directions))]

texttable <- data.table(Species = "Jack pine", 
                        Direction = directions,
                        x = 1989, y = Inf, 
                        texts = letters[1:length(directions)])
# texttable <- rbind(texttable, 
#                    data.table(Species = c("Jack pine", "Trembling aspen",
#                                           "Black spruce", "Other species"),
#                               Direction = "changewithIntraH",
#                               x = 1989, y = Inf, texts = " "))
texttable[,':='(Species = factor(Species, 
                                 levels = fourspecies),
                Direction = factor(Direction, 
                                   levels = directions))]
controPoints <- data.table(expand.grid(Direction = factor(c("mainTrend", "changewithInterH"),
                                                          levels = directions),
                                       PredictedABGR = c(0.026, 0.134, 0.687)))
controPoints[,':='(Species = factor("Jack pine", levels = fourspecies),
                           Year = 1995)]
controPoints <- controPoints[PredictedABGR != 0.026 | Direction != "mainTrend",]

figure <- ggplot(data = allFigureData[Direction != "mainTrend"], aes(x = Year, y = PredictedABGR))+
  facet_grid(Direction~Species, scale = "free_y",  drop = FALSE)+
  geom_point(data = controPoints, aes(x = Year, y = PredictedABGR), col = "white")+
  geom_line(aes(group = CompetitionIntensity, col = CompetitionIntensity), 
            size = 1)+
  geom_ribbon(data = allFigureData[Direction == "mainTrend", ], 
              aes (x = Year, ymin = PredictedABGR_Lower, 
                   ymax = PredictedABGR_Upper),
            fill = "gray", alpha = 0.5)+
  geom_line(data = allFigureData[Direction == "mainTrend", ], 
            aes (x = Year, y = PredictedABGR),
            col = "black", size = 1)+
  # geom_line(aes(group = CompetitionIntensity, y = PredictedABGR_Lower, col = CompetitionIntensity), linetype = 2)+
  # geom_line(aes(group = CompetitionIntensity, y = PredictedABGR_Upper, col = CompetitionIntensity), linetype = 2)+
  scale_y_log10(name = expression(paste("Aboveground biomass growth rate (Kg ", year^{-1}, ")")),
                breaks = c(0.001, 0.005, 0.026, 0.134, 0.687, 3.516, 18))+
  scale_x_continuous(name = "Year", breaks = seq(1990, 2010, by = 5))+
  scale_color_manual(name = "Competition intensity", values = c("blue", "black", "magenta"),
                     labels = c("Strong", "Medium strong", "Weak"))+
  guides(col = guide_legend(title.position = "top",
                            direction = "horizontal"))+
  geom_segment(data = segmenttable, aes(x = x, xend = xend, y = y, yend = yend), size = 1)+
  geom_text(data = texttable, aes(x = x, y = y, label = texts), vjust = 1, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # panel.margin = unit(0, "line"),
        # axis.line.x = element_line(size = 1, colour = "black"),
        # axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour = "white"),
        strip.text.x = element_text(size = 15),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.20, 0.85),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))


dev(4)
clearPlot()
ggsave(file = file.path(workPath, "TablesFigures", "Figure 3. temporal trends by competition.png"),
       figure,  width = 12, height = 8)






