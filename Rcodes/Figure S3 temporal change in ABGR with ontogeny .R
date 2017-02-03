rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
load(file.path(workPath, "data", selectionMethod, "BestYearModels.RData"))
##### for overall temporal trends and its dependency on DBH and RBI
analysesData[,SA:=IniFA+2.5]
##### for overall temporal trends and its dependency on DBH and RBI
for(i in 1:length(allHFixedCoeff)){
  indicoeff <- allHFixedCoeff[[i]]
  indicoeff[,Species:=names(allHFixedCoeff[i])]
  if(i == 1){
    allHcoeff <- indicoeff
  } else {
    allHcoeff <- rbind(allHcoeff, indicoeff)
  }
}
allHcoeff[, variable:= unlist(lapply(lapply(lapply((rn), function(x) unlist(strsplit(x, ":"))), function(y) sort(y)),
                                     function(z) paste(z, collapse = ":")))]
allHcoeff[variable == "logDBHctd:Yearctd", Direction:="changewithDBH"]
allHcoeff[variable == "logSActd:Yearctd", Direction:="changewithSA"]

output <- data.table(Species = character(), Direction = character(),
                     Year = numeric(),
                     CompetitionIntensity = numeric(), PredictedABGR = numeric(),
                     PredictedABGR_Lower = numeric(), PredictedABGR_Upper = numeric(),
                     overallSignificant = numeric())
for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies,]
  speciesallHcoeff <- allHcoeff[Species == indispecies, ]
  allHbestFormula <- allHbestFormulas[[indispecies]]
  theallHmodel <- allHbestModels[[indispecies]]
  # change with H
  if(nrow(speciesallHcoeff[Direction == "changewithDBH",]) == 1){
    changewithH <- data.table(expand.grid(Species = indispecies,
                                          Direction = "changewithDBH",
                                          Year = seq(min(speciesData$Year), 
                                                     max(speciesData$Year),
                                                     length = 100), 
                                          DBH = exp(seq(log(min(speciesData$IniDBH)),
                                                      log(max(speciesData$IniDBH)),
                                                      length = 100)),
                                          stringsAsFactors = FALSE))
    changewithH[,':='(Yearctd = Year-mean(speciesData$Year),
                      logDBHctd = log(DBH)-mean(log(speciesData$IniDBH)),
                      logHctd = 0,
                      logSActd = 0)]
    changewithH[,CompetitionIntensity:=as.numeric(as.factor(DBH))]
    fittedvalues <- predict(theallHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
    changewithH$PredictedABGR <- exp(fittedvalues$fit)
    changewithH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    
    output <- rbind(output, changewithH[,.(Species, Direction, Year, 
                                           CompetitionIntensity,
                                           PredictedABGR, PredictedABGR_Lower,
                                           PredictedABGR_Upper,
                                           overallSignificant = 0)])
    rm(fittedvalues, changewithH)
  }
  if(nrow(speciesallHcoeff[Direction == "changewithSA",]) == 1){
    changewithH <- data.table(expand.grid(Species = indispecies,
                                          Direction = "changewithSA",
                                          Year = seq(min(speciesData$Year), 
                                                     max(speciesData$Year),
                                                     length = 100), 
                                          SA = exp(seq(log(min(speciesData$SA)),
                                                        log(max(speciesData$SA)),
                                                        length = 100)),
                                          stringsAsFactors = FALSE))
    changewithH[,':='(Yearctd = Year-mean(speciesData$Year),
                      logDBHctd = 0,
                      logHctd = 0,
                      logSActd = log(SA)-mean(log(speciesData$SA)))]
    changewithH[,CompetitionIntensity:=as.numeric(as.factor(SA))]
    fittedvalues <- predict(theallHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
    changewithH$PredictedABGR <- exp(fittedvalues$fit)
    changewithH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
    changewithH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
    
    output <- rbind(output, changewithH[,.(Species, Direction, Year, 
                                           CompetitionIntensity,
                                           PredictedABGR, PredictedABGR_Lower,
                                           PredictedABGR_Upper,
                                           overallSignificant = 0)])
    rm(fittedvalues, changewithH)
  }
}

output[Species == "Minor species", Species:="Minor species group"]
output[Species == "All species", Species:="All trees"]
studySpecies <- c("All trees", studySpecies[2:4], "Minor species group")

output[,':='(Species = factor(Species, studySpecies),
             Direction = factor(Direction, 
                                levels = c("changewithDBH", "changewithSA")))]


segmenttable <- data.table(Species = "All trees", 
                           Direction = unique(output$Direction),
                           x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

segmenttable2 <- data.table(expand.grid(Species = factor(studySpecies),
                                        Direction = unique(output$Direction)))
segmenttable2[,':='(x = -Inf, xend = Inf, y = -Inf, yend = -Inf)]
segmenttable <- rbind(segmenttable, segmenttable2)



segmenttable[,':='(Species = factor(Species, 
                                    levels = studySpecies),
                   Direction = factor(Direction, 
                                      levels = c("changewithDBH", "changewithSA")))]

texttable <- data.table(Species = "All trees", 
                        Direction = unique(output$Direction),
                        x = 1989, y = Inf)
texttable[,':='(Species = factor(Species, 
                                 levels = studySpecies),
                Direction = factor(Direction, 
                                   levels = c("changewithDBH", "changewithSA")))]
texttable[,texts := letters[as.numeric(Direction)]]

figure <- ggplot(data = output, aes(x = Year, y = PredictedABGR))+
  facet_grid(Direction~Species, scale = "free_y",  drop = TRUE)+
  geom_line(aes(group = CompetitionIntensity, col = CompetitionIntensity), 
            size = 1)+
  scale_y_continuous(name = expression(paste("Tree-level aboveground biomass growth rate (Kg ", year^{-1}, ")")))+
  scale_x_continuous(name = "Year", breaks = seq(1990, 2010, by = 5))+
  scale_color_continuous(name = "DBH / SA",
                         low = "blue", high = "red", breaks = c(3, 100),
                         labels = c("\nSmall / young", "Big / old\n"),
                         guide = guide_colorbar(reverse = TRUE, ticks = FALSE))+
  geom_segment(data = segmenttable, aes(x = x, xend = xend, y = y, yend = yend), size = 1)+
  geom_text(data = texttable, aes(x = x, y = y, label = texts), vjust = 1, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_text(size = 15, face = "italic"),
        # legend.direction = "horizontal",
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.85, 0.35),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))
workPath <- "~/GitHub/Climate_Growth"
ggsave(file = file.path(workPath, "TablesFigures", 
                        paste("Figure S3. temporal trends with DBHSA", selectionMethod, ".png")),
       figure,  width = 11, height = 9.5)
