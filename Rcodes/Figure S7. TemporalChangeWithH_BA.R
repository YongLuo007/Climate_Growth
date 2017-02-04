rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
load(file.path(workPath, "data",selectionMethod,
               "BestYearModels_BA.RData"))
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
allHcoeff[variable == "logHctd:Yearctd", Direction:="changewithH"]

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
  mainTrends <- data.table(Species = indispecies, 
                           Direction = "mainTrend",
                           Year = seq(min(speciesData$Year), 
                                      max(speciesData$Year),
                                      length = 100))
  mainTrends[,':='(Yearctd = Year - mean(speciesData$Year),
                   logDBHctd = 0,
                   logSActd = 0,
                   logHctd = 0)]
  fittedvalues <- predict(theallHmodel, newdata = mainTrends, level = 0, se.fit = TRUE)
  treeToPlotConvertor <- (length(unique(speciesData$uniTreeID))/length(unique(speciesData$PlotID)))
  mainTrends$PredictedABGR <- treeToPlotConvertor*exp(fittedvalues$fit)
  mainTrends$PredictedABGR_Upper <- treeToPlotConvertor*exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  mainTrends$PredictedABGR_Lower <- treeToPlotConvertor*exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  if(nrow(speciesallHcoeff[rn == "Yearctd",]) == 1){
    output <- rbind(output, mainTrends[,.(Species, Direction, Year, 
                                          CompetitionIntensity = 0,
                                          PredictedABGR, PredictedABGR_Lower,
                                          PredictedABGR_Upper,
                                          overallSignificant = 1)])
    
  } else {
    output <- rbind(output, mainTrends[,.(Species, Direction, Year, 
                                          CompetitionIntensity = 0,
                                          PredictedABGR, PredictedABGR_Lower,
                                          PredictedABGR_Upper,
                                          overallSignificant = 2)])
  }
  rm(fittedvalues, mainTrends)
  
  # change with H
  if(nrow(speciesallHcoeff[Direction == "changewithH",]) == 1){
    changewithH <- data.table(expand.grid(Species = indispecies,
                                          Direction = "changewithH",
                                          Year = seq(min(speciesData$Year), 
                                                     max(speciesData$Year),
                                                     length = 100), 
                                          H = exp(seq(log(min(speciesData$H)),
                                                      log(max(speciesData$H)),
                                                      length = 100)),
                                          stringsAsFactors = FALSE))
    changewithH[,':='(Yearctd = Year-mean(speciesData$Year),
                      logDBHctd = 0,
                      logHctd = log(H)-mean(log(speciesData$H)),
                      logSActd = 0)]
    changewithH[,CompetitionIntensity:=as.numeric(as.factor(H))]
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

output[Species == "All species", Species:="All trees"]
output[Species == "Minor species", Species:="Minor species group"]
studySpecies <- c("All trees", studySpecies[2:4], "Minor species group")
output[,':='(Species = factor(Species, studySpecies),
             Direction = factor(Direction, 
                                levels = c("mainTrend", "changewithH")))]


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
                                      levels = c("mainTrend", "changewithH")))]

texttable <- data.table(Species = "All trees", 
                        Direction = unique(output$Direction),
                        x = 1989, y = Inf)
texttable[,':='(Species = factor(Species, 
                                 levels = studySpecies),
                Direction = factor(Direction, 
                                   levels = c("mainTrend", "changewithH")))]
texttable[,texts := letters[as.numeric(Direction)]]
newlabels1 <- list("All trees" = "All trees",
                   "Jack pine" = "Jack pine",
                   "Trembling aspen" = "Trembling aspen",
                   "Black spruce" = "Black spruce",
                   "Minor species group" = "Minor species group")
newlabels2 <- list("mainTrend" = expression(atop("Plot-level basal area growth rate",
                                                 paste("(",cm^{2}, " ",  year^{-1}, " per plot)"))),
                   "changewithH" = expression(atop("Tree-level basal area growth rate",
                                                   paste("(", cm^{2}, " ", year^{-1}, " per tree)"))))


figure_labeller <- function(variable,value){
  if(variable == "Species"){
    return(newlabels1[value])
  } else if (variable == "Direction"){
    return(newlabels2[value])
  } 
}



figure <- ggplot(data = output[Direction != "mainTrend"], aes(x = Year, y = PredictedABGR))+
  facet_grid(Direction~Species,  drop = TRUE,
             switch = "y",
             labeller = figure_labeller)+
  geom_line(aes(group = CompetitionIntensity, col = CompetitionIntensity), 
            size = 1)+
  geom_ribbon(data = output[Direction == "mainTrend", ], 
              aes(x = Year, ymin = PredictedABGR_Lower, ymax = PredictedABGR_Upper),
              col = "white", fill = "gray", show.legend = FALSE)+
  geom_line(data = output[Direction == "mainTrend", ], 
            aes(x = Year, y = PredictedABGR,
                linetype = as.factor(overallSignificant)),
            col = "black", size = 1, show.legend = FALSE)+
  # scale_y_continuous(name = expression(paste("Aboveground biomass growth rate (Kg ", year^{-1}, ")")))+
  scale_x_continuous(name = "Year", breaks = seq(1990, 2010, by = 5))+
  scale_color_continuous(name = "Competition \nintensity",
                         low = "blue", high = "red", breaks = c(3, 100),
                         labels = c("\nWeak", "Strong\n"),
                         guide = guide_colorbar(reverse = TRUE, ticks = FALSE))+
  geom_segment(data = segmenttable, aes(x = x, xend = xend, y = y, yend = yend), size = 1)+
  geom_text(data = texttable, aes(x = x, y = y, label = texts), vjust = 1, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_text(size = 15, face = "italic"),
        strip.text.y = element_text(size = 15),
        # legend.direction = "horizontal",
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.85, 0.35),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))
workPath <- "~/GitHub/Climate_Growth"
ggsave(file = file.path(workPath, "TablesFigures", 
                        paste("Figure S7. temporal trends in BA with H.png")),
       figure,  width = 11, height = 9.5)
