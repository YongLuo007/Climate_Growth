rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "Results", "bestYearModel.RData"))
rm(speciesData, indispecies)
##### for overall temporal trends and its dependency on DBH and RBI
studySpecies <- c("All", "JP", "TA", "BS", "Other")
# competition effect change with year
tempallCoeff <- lapply(allbestmodels, function(x){data.table(summary(x)$tTable, keep.rownames = TRUE)})
for(i in 1:length(tempallCoeff)){
  coeffIndi <- cbind(data.table(Species = rep(names(tempallCoeff)[i], nrow(tempallCoeff[[i]]))), 
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

# figure 2 competition effect change with dominance
CompetitionWithYearTable <- data.table(Species = character(), 
                                       Year = numeric(),
                                       Value = numeric(), SE = numeric(), MainEffect = numeric())
for(indispecies in studySpecies){
  set.seed(1)
  if(indispecies == "All"){
    speciesData <- data.table::copy(analysesData)
  } else if(indispecies == "Other"){
    speciesData <- analysesData[!(Species %in% c("JP", "TA", "BS")),]
  } else {
    speciesData <- analysesData[Species == indispecies,]
  }
  YearRange <- range(speciesData$Year)
  CompetitionWithYearOutput_indispecies <- data.table(Year = seq(YearRange[1], YearRange[2], length = 20))[
    , Yearctd:=Year-mean(speciesData$Year)][,.(Species = indispecies, Year, RBI = 0, Yearctd)]
  CompetitionWithYearDomOutput_indispecies <- data.table(expand.grid(Year = seq(YearRange[1], YearRange[2], length = 100),
                                     RBI = c(seq(0, 100, length = 10))))[
                                       , ':='(Yearctd = Year - mean(speciesData$Year),
                                              RBIctd = RBI - mean(speciesData$RBI))][
                                                ,.(Species = indispecies, Year, RBI, Yearctd, RBIctd)]
  coeffTable <- allCoeff[Species == indispecies,]
  ## get intra h and interh main effect
  if(nrow(coeffTable[rn == "logIntraHctd",])>0){
    MainIntraHeffect <- rnorm(10000,
                          mean = coeffTable[rn == "logIntraHctd",]$Value,
                          sd = 100*coeffTable[rn == "logIntraHctd",]$Std.Error)
  } else {
    MainIntraHeffect <- 0
  }
  if(nrow(coeffTable[rn == "logInterHctd",])>0){
    MainInterHeffect <- rnorm(10000,
                          mean = coeffTable[rn == "logInterHctd",]$Value,
                          sd = 100*coeffTable[rn == "logInterHctd",]$Std.Error)
  } else {
    MainInterHeffect <- 0
  }
  # intraH and year, interH and year interaction
  if(nrow(coeffTable[rn == "logIntraHctd:Yearctd",])>0){
    IntraHwithYear <- rnorm(10000,
                            mean = coeffTable[rn == "logIntraHctd:Yearctd",]$Value,
                            sd = 100*coeffTable[rn == "logIntraHctd:Yearctd",]$Std.Error)
  } else {
    IntraHwithYear <- 0
  }
  if(nrow(coeffTable[rn == "logInterHctd:Yearctd",])>0){
    InterHwithYear <- rnorm(10000,
                            mean = coeffTable[rn == "logInterHctd:Yearctd",]$Value,
                            sd = 100*coeffTable[rn == "logInterHctd:Yearctd",]$Std.Error)
  } else {
    InterHwithYear <- 0
  }
  # intraH and RBI, interH and RBI interaction
  if(nrow(coeffTable[rn == "logIntraHctd:RBIctd",])>0){
    IntraHwithRBI <- rnorm(10000,
                            mean = coeffTable[rn == "logIntraHctd:RBIctd",]$Value,
                            sd = 100*coeffTable[rn == "logIntraHctd:RBIctd",]$Std.Error)
  } else {
    IntraHwithRBI <- 0
  }
  if(nrow(coeffTable[rn == "logInterHctd:RBIctd",])>0){
    InterHwithRBI <- rnorm(10000,
                           mean = coeffTable[rn == "logInterHctd:RBIctd",]$Value,
                           sd = 100*coeffTable[rn == "logInterHctd:RBIctd",]$Std.Error)
  } else {
    InterHwithRBI <- 0
  }
  # tree way interactions
  if(nrow(coeffTable[rn == "logIntraHctd:RBIctd:Yearctd",])>0){
    IntraHwithYearRBI <- rnorm(10000,
                            mean = coeffTable[rn == "logIntraHctd:RBIctd:Yearctd",]$Value,
                            sd = 100*coeffTable[rn == "logIntraHctd:RBIctd:Yearctd",]$Std.Error)
  } else {
    IntraHwithYearRBI <- 0
  }
  if(nrow(coeffTable[rn == "logInterHctd:RBIctd:Yearctd",])>0){
    InterHwithYearRBI <- rnorm(10000,
                            mean = coeffTable[rn == "logInterHctd:RBIctd:Yearctd",]$Value,
                            sd = 100*coeffTable[rn == "logInterHctd:RBIctd:Yearctd",]$Std.Error)
  } else {
    InterHwithYearRBI <- 0
  }
  for(i in 1:nrow(CompetitionWithYearOutput_indispecies)){
    CompetitionWithYearOutput_indispecies$IntraH_Value[i] <- 
      mean(CompetitionWithYearOutput_indispecies$Yearctd[i]*
             IntraHwithYear+MainIntraHeffect)
    CompetitionWithYearOutput_indispecies$IntraH_SE[i] <- 
      sd(CompetitionWithYearOutput_indispecies$Yearctd[i]*
           IntraHwithYear+MainIntraHeffect)/100
    CompetitionWithYearOutput_indispecies$InterH_Value[i] <- 
      mean(CompetitionWithYearOutput_indispecies$Yearctd[i]*
             InterHwithYear+MainInterHeffect)
    CompetitionWithYearOutput_indispecies$InterH_SE[i] <- 
      sd(CompetitionWithYearOutput_indispecies$Yearctd[i]*
           InterHwithYear+MainInterHeffect)/100
  }
  for(i in 1:nrow(CompetitionWithYearDomOutput_indispecies)){
    CompetitionWithYearDomOutput_indispecies$IntraH_Value[i] <- 
      mean(MainIntraHeffect+CompetitionWithYearDomOutput_indispecies$Yearctd[i]*
      IntraHwithYear+CompetitionWithYearDomOutput_indispecies$RBIctd[i]*
      IntraHwithRBI+CompetitionWithYearDomOutput_indispecies$Yearctd[i]*
      CompetitionWithYearDomOutput_indispecies$RBIctd[i]*IntraHwithYearRBI)
    CompetitionWithYearDomOutput_indispecies$IntraH_SE[i] <- 
      sd(MainIntraHeffect+CompetitionWithYearDomOutput_indispecies$Yearctd[i]*
             IntraHwithYear+CompetitionWithYearDomOutput_indispecies$RBIctd[i]*
             IntraHwithRBI+CompetitionWithYearDomOutput_indispecies$Yearctd[i]*
             CompetitionWithYearDomOutput_indispecies$RBIctd[i]*IntraHwithYearRBI)
    CompetitionWithYearDomOutput_indispecies$InterH_Value[i] <- 
      mean(MainInterHeffect+CompetitionWithYearDomOutput_indispecies$Yearctd[i]*
             InterHwithYear+CompetitionWithYearDomOutput_indispecies$RBIctd[i]*
             InterHwithRBI+CompetitionWithYearDomOutput_indispecies$Yearctd[i]*
             CompetitionWithYearDomOutput_indispecies$RBIctd[i]*InterHwithYearRBI)
    CompetitionWithYearDomOutput_indispecies$InterH_SE[i] <- 
      sd(MainInterHeffect+CompetitionWithYearDomOutput_indispecies$Yearctd[i]*
             InterHwithYear+CompetitionWithYearDomOutput_indispecies$RBIctd[i]*
             InterHwithRBI+CompetitionWithYearDomOutput_indispecies$Yearctd[i]*
             CompetitionWithYearDomOutput_indispecies$RBIctd[i]*InterHwithYearRBI)
    
  }
  CompetitionWithYearOutput_indispecies <- 
    rbind(CompetitionWithYearOutput_indispecies[,Main:=0],
    data.table(Species = indispecies, Year = mean(speciesData$Year), RBI = 0, Yearctd = 0, 
               IntraH_Value = coeffTable[rn == "logIntraHctd",]$Value,
               IntraH_SE = coeffTable[rn == "logIntraHctd",]$Std.Error,
               InterH_Value = coeffTable[rn == "logInterHctd",]$Value,
               InterH_SE = coeffTable[rn == "logInterHctd",]$Std.Error,
               Main = 1))
  if(indispecies == "All"){
    CompetitionWithYearOutput <- CompetitionWithYearOutput_indispecies
    
    CompetitionWithYearDomOutput <- CompetitionWithYearDomOutput_indispecies
  } else {
    CompetitionWithYearOutput <- rbind(CompetitionWithYearOutput,
                                       CompetitionWithYearOutput_indispecies)
    CompetitionWithYearDomOutput <- rbind(CompetitionWithYearDomOutput,
                                          CompetitionWithYearDomOutput_indispecies)
  }
}
CompetitionWithYearOutput <- 
  rbind(CompetitionWithYearOutput[,.(Direction = "Overall", Competition = "IntraH",
                                     Species, Year, RBI, Value = IntraH_Value, SE = IntraH_SE, Main)],
        CompetitionWithYearOutput[,.(Direction = "Overall", Competition = "InterH",
                                     Species, Year, RBI, Value = InterH_Value, SE = InterH_SE, Main)])
CompetitionWithYearDomOutput <- 
  rbind(CompetitionWithYearDomOutput[,.(Direction = "IntraH", Competition = "IntraH",
                                        Species, Year, RBI, Value = IntraH_Value, SE = IntraH_SE, Main = 0)],
        CompetitionWithYearDomOutput[,.(Direction = "InterH", Competition = "InterH",
                                        Species, Year, RBI, Value = InterH_Value, SE = InterH_SE, Main = 0)])

figureData <- rbind(CompetitionWithYearDomOutput[, SE:=0], CompetitionWithYearOutput[Main == 1,])
figureData[,':='(Direction = factor(Direction, levels = c("Overall", "IntraH", "InterH"),
                                    labels = c("a", "b", "c")),
                 Species = factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                                  labels = c("All species", "Jack pine", "Trembling aspen", 
                                             "Black spruce", "Other species")),
                 Competition = factor(Competition, levels = c("IntraH", "InterH")))]
# values <- figureData[,.(minV = min(Value-1.99*SE), maxV = max(Value+1.99*SE)), by = Species]
# values <- reshape(values, varying = c("minV", "maxV"), v.names = "Value", direction = "long") %>%
#   data.table
# dummyData <- data.table(Direction = factor(c("a", "b", "c")), k = 1, key = "k")[
#   setkey(values[,.(Species, Value, k = 1)], k), allow.cartesian = TRUE][
#     , ':='(k = NULL, Competition = 1.5, Year = 1986)]

equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    round(seq(min(x)+d, max(x)-d, length=n), 2)
  }
}  

figureData1 <- figureData[Direction %in% c("b", "c"),][
  !(Species %in% c("Jack pine", "Trembling aspen", "Black spruce", "Other species") &
      Direction == "c"),][!(Species == "Other species" & Direction == "b"),]
yaxisLineData <- data.table(expand.grid(Direction = c("b", "c"), 
                                        Species = c("All species", "Jack pine", "Trembling aspen", 
                                                    "Black spruce", "Other species")))[
                                                      ,':='(Year = -Inf, Yearend = -Inf, 
                                                      Value = Inf, Valueend = -Inf)]
mainEffect <- figureData[Main == 1, ][, ':='(Direction = factor(Competition, 
                                                                levels = c("IntraH", "InterH"),
                                                                labels = c("b", "c")))]
Figure3 <- ggplot(data = figureData1[Direction %in% c("b", "c"),], aes(x = Year, y = Value))+
  facet_grid(Species~Direction)+
  geom_line(aes(col = RBI, group = RBI), size = 1)+
  geom_point(data = mainEffect, aes(x = Year, y = Value), col = "blue", size = 2)+
  geom_errorbar(data = mainEffect, aes(ymin = Value-1.98*SE, ymax = Value+1.98*SE), 
                width = 1, col = "blue", size = 1)+
  scale_colour_continuous(name = "RBI", low = "#FF0000", high = "#00FF00", breaks = c(0, 50, 100))+
  facet_grid(Species~Direction, scales = "free_y")+
  geom_segment(data = yaxisLineData,
               aes(x = Year, xend = Yearend, y = Value, yend = Valueend), size = 1.5)+
  geom_text(data = data.frame(Direction = c("b", "c"), labels = c("a", "b"), 
                              Species = "All species",
                              Value = -0.05, Year=1986),
            aes(x = Year, y = Value, label = labels), size = 10)+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  scale_y_continuous(name = "Effect of competition on tree growth", 
                     breaks = equal_breaks(n=4, s=0.05))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.background = element_rect(colour = "white"),
        strip.text.y = element_text(size = 16),
        strip.text.x = element_blank(),
        legend.position = c(0.6, 0.6),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

dev(4)
clearPlot()

ggsave(file = file.path(workPath, "TablesFigures", "Figure 3_temporal trends in competition.png"),
       Figure3,  width = 10, height = 10)






