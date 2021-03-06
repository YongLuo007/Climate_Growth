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
allcoeff[variable == "logSActd:Yearctd", Direction:="changewithSA"]
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


output[,':='(Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen",
                                                  "Black spruce", "Other species")))]
plotFigure <- "OntogenyOnly"
if(plotFigure == "CompetitionOnly"){
  directions <- c("mainTrend", "changewithIntraH", "changewithInterH")
} else if (plotFigure == "OntogenyOnly") {
  directions <- c("changewithDBH", "changewithSA")
} else {
  directions <- c("mainTrend", "changewithIntraH", "changewithInterH", "changewithDBH", "changewithSA")
}
allFigureData <- data.table::copy(output)[Direction %in% directions, ]
allFigureData[, ':='(Direction = factor(Direction, levels = directions),
                     CompetitionIntensity = factor(CompetitionIntensity, 
                                                   levels = c("Strong", "Medium", "Weak")))]

segmenttable <- data.table(Species = "Jack pine", 
                           Direction = directions,
                           x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
segmenttable2 <- data.table(expand.grid(Species = factor(c("Jack pine", "Trembling aspen",
                                                           "Black spruce", "Other species")),
                                        Direction = directions))
segmenttable <- rbind(segmenttable, segmenttable2[,':='(x = -Inf, xend = Inf, y = -Inf, yend = -Inf)])
segmenttable[,':='(Species = factor(Species, 
                                    levels = c("Jack pine", "Trembling aspen",
                                               "Black spruce", "Other species")),
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
                                 levels = studySpecies),
                Direction = factor(Direction, 
                                   levels = directions))]
legendTableb <- texttable[Direction == "changewithSA",.(Species, Direction, x = 1998,
                                                        y = c(3.2, 2.9, 2.6),
                                                        texts = c("Maximum", "Medium",
                                                                  "Minimum"))]
legendTablebLine <- legendTableb[,.(Species, Direction, x = 1995, xend = 1997, y, yend = y)]
legendTablebTitle <- texttable[Direction == "changewithSA", .(Species, Direction, x = 1995,
                                                              y = 3.5, texts = "Stand age (SA)")]
legendTablebRect <- texttable[Direction == "changewithSA", .(Species, Direction, Year = 1994, xmax = 2009,
                                                             PredictedABGR = 2.38, ymax = 3.65)]
figure <- ggplot(data = allFigureData[Direction != "mainTrend"], aes(x = Year, y = PredictedABGR))+
  facet_grid(Direction~Species, scales = "free_y")+
  geom_line(aes(group = CompetitionIntensity, col = CompetitionIntensity), 
            size = 1)+
  # geom_line(aes(group = CompetitionIntensity, y = PredictedABGR_Lower, col = CompetitionIntensity), linetype = 2)+
  # geom_line(aes(group = CompetitionIntensity, y = PredictedABGR_Upper, col = CompetitionIntensity), linetype = 2)+
  scale_y_continuous(name = expression(paste("Aboveground biomass growth rate (Kg ", year^{-1}, ")")))+
  scale_x_continuous(name = "Year", breaks = seq(1990, 2010, by = 5))+
  scale_color_manual(name = "Diameter at breast height (DBH)", values = c("blue", "black", "magenta"),
                     labels = c("Maximum", "Medium", "Minimum"))+
  guides(col = guide_legend(title.position = "top"))+
  geom_segment(data = segmenttable, aes(x = x, xend = xend, y = y, yend = yend), size = 1)+
  geom_text(data = texttable, aes(x = x, y = y, label = texts), vjust = 1.5, size = 10)+
  geom_text(data = legendTableb, aes(x = x, y = y, label = texts), hjust = 0, size = 5)+
  geom_text(data = legendTablebTitle, aes(x = x, y = y, label = texts), hjust = 0, size = 6)+
  geom_rect(data = legendTablebRect, aes(xmin = Year, xmax = xmax, ymin = PredictedABGR, ymax = ymax),
            colour = "black", alpha = 0)+
  geom_segment(data = legendTablebLine, aes(x = x, xend = xend, y = y, yend = yend), size = 1,
               col = c("blue", "black", "magenta"))+
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
        legend.key = element_rect(colour = "white"),
        legend.position = c(0.20, 0.85),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

ggsave(file = file.path(workPath, "TablesFigures", "Figure S3. temporal trends by ontogeny.png"),
       figure,  width = 12, height = 8)






