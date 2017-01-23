rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data",
               "finalYearModels_AllCensus_PositiveGrowth.RData"))
##### for overall temporal trends and its dependency on DBH and RBI
for(i in 1:length(allHFixedCoeff)){
  indicoeff <- allHFixedCoeff[[i]]
  indicoeff[,Species:=names(allHFixedCoeff[i])]
  indicoeff2 <- indiHFixedCoeff[[i]]
  indicoeff2[,Species:=names(allHFixedCoeff[i])]
  if(i == 1){
    allHcoeff <- indicoeff
    indiHcoeff <- indicoeff2
  } else {
    allHcoeff <- rbind(allHcoeff, indicoeff)
    indiHcoeff <- rbind(indiHcoeff, indicoeff2)
  }
}
allHcoeff[, variable:= unlist(lapply(lapply(lapply((rn), function(x) unlist(strsplit(x, ":"))), function(y) sort(y)),
                                     function(z) paste(z, collapse = ":")))]
allHcoeff[variable == "logHctd:Yearctd", Direction:="changewithH"]
indiHcoeff[, variable:= unlist(lapply(lapply(lapply((rn), function(x) unlist(strsplit(x, ":"))), function(y) sort(y)),
                                      function(z) paste(z, collapse = ":")))]
indiHcoeff[variable == "logIntraHctd:Yearctd", Direction:="changewithIntraH"]
indiHcoeff[variable == "logInterHctd:Yearctd", Direction:="changewithInterH"]

output <- data.table(Model = character(), Species = character(), Direction = character(), Year = numeric(),
                     CompetitionIntensity = numeric(), PredictedABGR = numeric(),
                     PredictedABGR_Lower = numeric(), PredictedABGR_Upper = numeric())
for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies,]
  speciesallHcoeff <- allHcoeff[Species == indispecies, ]
  speciesindiHcoeff <- indiHcoeff[Species == indispecies,]
  allHbestFormula <- allHbestFormulas[[indispecies]]
  theallHmodel <- allHbestModels[[indispecies]]
  indiHbestFormu <- indiHbestFormulas[[indispecies]]
  theindiHmodel <- indiHbestModels[[indispecies]]
  if(nrow(speciesallHcoeff[rn == "Yearctd",]) == 1){
    mainTrends <- data.table(Model = "allH",
                             Species = indispecies, 
                             Direction = "mainTrend",
                             Year = seq(min(speciesData$Year), 
                                        max(speciesData$Year),
                                        length = 100))
    mainTrends[,':='(Yearctd = Year - mean(speciesData$Year),
                     logDBHctd = 0,
                     logSActd = 0,
                     logHctd = 0)]
    fittedvalues <- predict(theallHmodel, newdata = mainTrends, level = 0, se.fit = TRUE)
    mainTrends$PredictedABGR <- exp(fittedvalues$fit)
    mainTrends$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    mainTrends$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    output <- rbind(output, mainTrends[,.(Model, Species, Direction, Year, 
                                          CompetitionIntensity = 0,
                                          PredictedABGR, PredictedABGR_Lower,
                                          PredictedABGR_Upper)])
    rm(fittedvalues, mainTrends)
  }
  if(nrow(speciesindiHcoeff[rn == "Yearctd",]) == 1){
    mainTrends <- data.table(Model = "indiH",
                             Species = indispecies, 
                             Direction = "mainTrend",
                             Year = seq(min(speciesData$Year), 
                                        max(speciesData$Year),
                                        length = 100))
    mainTrends <- rbind(data.table::copy(mainTrends)[,Model := paste(Model, "_IntraH", sep = "")],
                        data.table::copy(mainTrends)[,Model := paste(Model, "_InterH", sep = "")])
    mainTrends[,':='(Yearctd = Year - mean(speciesData$Year),
                     logDBHctd = 0,
                     logSActd = 0,
                     logIntraHctd = 0, logInterHctd = 0)]
    fittedvalues <- predict(theindiHmodel, newdata = mainTrends, level = 0, se.fit = TRUE)
    mainTrends$PredictedABGR <- exp(fittedvalues$fit)
    mainTrends$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    mainTrends$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    output <- rbind(output, mainTrends[,.(Model, Species, Direction, Year, 
                                          CompetitionIntensity = 0,
                                          PredictedABGR, PredictedABGR_Lower,
                                          PredictedABGR_Upper)])
    rm(fittedvalues, mainTrends)
  }
  # change with H
  if(nrow(speciesallHcoeff[Direction == "changewithH",]) == 1){
    changewithH <- data.table(expand.grid(Species = indispecies,
                                          Direction = "changewithH",
                                          Year = seq(min(speciesData$Year), 
                                                     max(speciesData$Year),
                                                     length = 100), 
                                          H = exp(seq(log(min(speciesData$H+1)),
                                                      log(max(speciesData$H+1)),
                                                      length = 100))-1,
                                          stringsAsFactors = FALSE))
    changewithH[,':='(Yearctd = Year-mean(speciesData$Year),
                      logDBHctd = 0,
                      logHctd = log(H+1)-mean(log(speciesData$H+1)),
                      logSActd = 0)]
    changewithH[,CompetitionIntensity:=as.numeric(as.factor(H))]
    
    
    fittedvalues <- predict(theallHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
    changewithH$PredictedABGR <- exp(fittedvalues$fit)
    changewithH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    
    output <- rbind(output, changewithH[,.(Model = "allH",
                                           Species, Direction, Year, 
                                           CompetitionIntensity,
                                           PredictedABGR, PredictedABGR_Lower, PredictedABGR_Upper)])
    rm(fittedvalues, changewithH)
  }
  # chnage with intraH
  
  if(nrow(speciesindiHcoeff[Direction == "changewithIntraH",]) == 1){
    changewithH <- data.table(expand.grid(Species = indispecies,
                                          Direction = "changewithIntraH",
                                          Year = seq(min(speciesData$Year), 
                                                     max(speciesData$Year),
                                                     length = 100), 
                                          IntraH = exp(seq(log(min(speciesData$IntraH+1)),
                                                           log(max(speciesData$IntraH+1)),
                                                           length = 100))-1,
                                          stringsAsFactors = FALSE))
    changewithH[,':='(Yearctd = Year-mean(speciesData$Year),
                      logDBHctd = 0,
                      logIntraHctd = log(IntraH+1)-mean(log(speciesData$IntraH+1)),
                      logInterHctd = 0,
                      logSActd = 0)]
    changewithH[,CompetitionIntensity:=as.numeric(as.factor(logIntraHctd))]
    fittedvalues <- predict(theindiHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
    changewithH$PredictedABGR <- exp(fittedvalues$fit)
    changewithH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    
    output <- rbind(output, changewithH[,.(Model = "indiH_IntraH",
                                           Species, Direction, Year, 
                                           CompetitionIntensity,
                                           PredictedABGR, PredictedABGR_Lower, PredictedABGR_Upper)])
    rm(fittedvalues, changewithH)
  }
  if(nrow(speciesindiHcoeff[Direction == "changewithInterH",]) == 1){
    changewithH <- data.table(expand.grid(Species = indispecies,
                                          Direction = "changewithInterH",
                                          Year = seq(min(speciesData$Year), 
                                                     max(speciesData$Year),
                                                     length = 100), 
                                          InterH = exp(seq(log(min(speciesData$InterH+1)),
                                                           log(max(speciesData$InterH+1)),
                                                           length = 100))-1,
                                          stringsAsFactors = FALSE))
    changewithH[,':='(Yearctd = Year-mean(speciesData$Year),
                      logDBHctd = 0,
                      logIntraHctd = 0,
                      logInterHctd = log(InterH+1)-mean(log(speciesData$InterH+1)),
                      logSActd = 0)]
    changewithH[,CompetitionIntensity:=as.numeric(as.factor(logInterHctd))]
    fittedvalues <- predict(theindiHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
    changewithH$PredictedABGR <- exp(fittedvalues$fit)
    changewithH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    
    output <- rbind(output, changewithH[,.(Model = "indiH_InterH",
                                           Species, Direction, Year, 
                                           CompetitionIntensity,
                                           PredictedABGR, PredictedABGR_Lower,
                                           PredictedABGR_Upper)])
    rm(fittedvalues)
  }
}


output[,':='(Species = factor(Species, studySpecies),
             Direction = factor(Direction, 
                                levels = c("mainTrend", "changewithH", "changewithIntraH", "changewithInterH")))]


segmenttable <- data.table(Species = "All species", 
                           Direction = unique(output$Direction),
                           x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

segmenttable2 <- data.table(expand.grid(Species = factor(studySpecies),
                                        Direction = unique(output$Direction)))
segmenttable2[,':='(x = -Inf, xend = Inf, y = -Inf, yend = -Inf)]
segmenttable <- rbind(segmenttable, segmenttable2)



segmenttable[,':='(Species = factor(Species, 
                                    levels = studySpecies),
                   Direction = factor(Direction, 
                                      levels = c("mainTrend", "changewithH", "changewithIntraH", "changewithInterH")))]

texttable <- data.table(Species = "All species", 
                        Direction = unique(output$Direction),
                        x = 1989, y = Inf)[,texts := letters[1:length(Direction)]]
texttable[,':='(Species = factor(Species, 
                                 levels = studySpecies),
                Direction = factor(Direction, 
                                   levels = c("mainTrend", "changewithH", "changewithIntraH", "changewithInterH")))]
figure <- ggplot(data = output[Direction != "mainTrend"], aes(x = Year, y = PredictedABGR))+
  facet_grid(Direction~Species, scale = "free_y",  drop = TRUE)+
  geom_line(aes(group = CompetitionIntensity, col = CompetitionIntensity), 
            size = 1)+
  geom_line(data = output[Direction == "mainTrend" & Model == "indiH_IntraH", ], 
            aes(x = Year, y = PredictedABGR),
            col = "gray", size = 1, show.legend = FALSE)+
  geom_line(data = output[Direction == "mainTrend" & Model == "allH", ], 
            aes(x = Year, y = PredictedABGR),
            col = "black", size = 1, show.legend = FALSE)+
  scale_y_continuous(name = expression(paste("Aboveground biomass growth rate (Kg ", year^{-1}, ")")))+
  scale_x_continuous(name = "Year", breaks = seq(1990, 2010, by = 5))+
  scale_color_continuous(name = "Competition \nintensity\n",
                       low = "blue", high = "red", breaks = c(3, 100),
                       labels = c("Weak", "Strong\n"),
                       guide = guide_colorbar(reverse = TRUE))+
  geom_segment(data = segmenttable, aes(x = x, xend = xend, y = y, yend = yend), size = 1)+
  geom_text(data = texttable, aes(x = x, y = y, label = texts), vjust = 1, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour = "white"),
        strip.text.x = element_text(size = 15),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.90, 0.15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))
workPath <- "~/GitHub/Climate_Growth"
ggsave(file = file.path(workPath, "TablesFigures", "Figure 3. temporal trends by competition.png"),
       figure,  width = 11, height = 9.5)






