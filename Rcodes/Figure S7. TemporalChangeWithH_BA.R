rm(list = ls())
# produce figure s7
library(data.table); 
library(ggplot2, 
        lib.loc = "~/GitHub/Climate_Growth/RequiredRPackages")
library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data",selectionMethod,
               "FullYearModels_BA.RData"))
##### for overall temporal trends and its dependency on DBH and RBI
for(i in 1:length(fixedCoeffAll)){
  indicoeff <- fixedCoeffAll[[i]]
  indicoeff[,Species:=names(fixedCoeffAll[i])]
  if(i == 1){
    allHcoeff <- indicoeff
  } else {
    allHcoeff <- rbind(allHcoeff, indicoeff)
  }
}
useLogQunatile95 <- TRUE
allHcoeff[, variable:= unlist(lapply(lapply(lapply((rn), function(x) unlist(strsplit(x, ":"))), function(y) sort(y)),
                                     function(z) paste(z, collapse = ":")))]
allHcoeff[variable == "logHctd:Yearctd", Direction:="changewithH"]

output <- data.table(Species = character(), Direction = character(),
                     Year = numeric(),
                     CompetitionIntensity = numeric(), PredictedBAGR = numeric(),
                     PredictedBAGR_Lower = numeric(), PredictedBAGR_Upper = numeric(),
                     overallSignificant = numeric())
for(indispecies in studySpecies){
  speciesData <- analysesDataAll[Species_Group == indispecies,]
  minABGR <- round(abs(min(speciesData$BAGR)), 3)+0.01
  speciesallHcoeff <- allHcoeff[Species == indispecies, ]
  theallHmodel <- allfullModelsAll[[indispecies]]
  mainTrends <- data.table(Species = indispecies, 
                           Direction = "mainTrend",
                           Year = seq(min(speciesData$IniYear), 
                                      max(speciesData$FinYear),
                                      length = 100))
  mainTrends[,':='(Yearctd = Year - mean(speciesData$MidYear),
                   logDBHctd = 0,
                   logSActd = 0,
                   logHctd = 0)]
  fittedvalues <- predict(theallHmodel, newdata = mainTrends, level = 0, se.fit = TRUE)
  treeToPlotConvertor <- (length(unique(speciesData$uniTreeID))/length(unique(analysesDataAll$PlotID)))
  mainTrends$PredictedBAGR <- treeToPlotConvertor*(exp(fittedvalues$fit)-minABGR)
  mainTrends$PredictedBAGR_Upper <- treeToPlotConvertor*(exp(fittedvalues$fit+1.98*fittedvalues$se.fit)-minABGR)
  mainTrends$PredictedBAGR_Lower <- treeToPlotConvertor*(exp(fittedvalues$fit-1.98*fittedvalues$se.fit)-minABGR)
  if(speciesallHcoeff[rn == "Yearctd",]$`p-value` < 0.05){
    output <- rbind(output, mainTrends[,.(Species, Direction, Year, 
                                          CompetitionIntensity = 0,
                                          PredictedBAGR, PredictedBAGR_Lower,
                                          PredictedBAGR_Upper,
                                          overallSignificant = 1)])
    
  } else {
    output <- rbind(output, mainTrends[,.(Species, Direction, Year, 
                                          CompetitionIntensity = 0,
                                          PredictedBAGR, PredictedBAGR_Lower,
                                          PredictedBAGR_Upper,
                                          overallSignificant = 2)])
  }
  rm(fittedvalues, mainTrends)
  
  # change with H
  
  if(nrow(speciesallHcoeff[Direction == "changewithH",]) == 1){
    if(useLogQunatile95){
      H95quantilelower <- exp(as.numeric(quantile(log(speciesData$H), probs = 0.025)))
      H95quantileupper <- exp(as.numeric(quantile(log(speciesData$H), probs = 0.975)))
      changewithH <- data.table(expand.grid(Species = indispecies,
                                            Direction = "changewithH",
                                            Year = seq(min(speciesData$IniYear), 
                                                       max(speciesData$FinYear),
                                                       length = 100), 
                                            H = exp(seq(log(H95quantilelower),
                                                        log(H95quantileupper),
                                                        length = 100)),
                                            stringsAsFactors = FALSE))
    } else {
      changewithH <- data.table(expand.grid(Species = indispecies,
                                            Direction = "changewithH",
                                            Year = seq(min(speciesData$IniYear), 
                                                       max(speciesData$FinYear),
                                                       length = 100), 
                                            H = exp(seq(log(min(speciesData$H)),
                                                        log(max(speciesData$H)),
                                                        length = 100)),
                                            stringsAsFactors = FALSE))
    }
    
    changewithH[,':='(Yearctd = Year-mean(speciesData$MidYear),
                      logDBHctd = 0,
                      logHctd = log(H)-mean(log(speciesData$H)),
                      logSActd = 0)]
    changewithH[,CompetitionIntensity:=as.numeric(as.factor(H))]
    fittedvalues <- predict(theallHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
    changewithH$PredictedBAGR <- exp(fittedvalues$fit)-minABGR
    changewithH$PredictedBAGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)-minABGR
    changewithH$PredictedBAGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)-minABGR
    
    output <- rbind(output, changewithH[,.(Species, Direction, Year, 
                                           CompetitionIntensity,
                                           PredictedBAGR, PredictedBAGR_Lower,
                                           PredictedBAGR_Upper,
                                           overallSignificant = 0)])
    rm(fittedvalues, changewithH)
  }
}

output[Species == "Minor species", Species:="Minor species group"]
output[Species == "All species", Species:="All trees"]

studySpecies <- c("All trees", studySpecies[2:4], "Minor species group")
output[,':='(Species = factor(Species, levels = studySpecies),
             Direction = factor(Direction, 
                                levels = c("changewithH", "mainTrend")))]

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
                                      levels = c("changewithH", "mainTrend")))]
texttable <- data.table(Species = "All trees", 
                        Direction = unique(output$Direction),
                        x = 1989, y = Inf)
texttable[,':='(Species = factor(Species, 
                                 levels = studySpecies),
                Direction = factor(Direction, 
                                   levels = c("changewithH", "mainTrend")))]
texttable[,texts := letters[as.numeric(Direction)]]

newlabels1 <- list("All trees" = "All trees",
                   "Jack pine" = "Jack pine",
                   "Trembling aspen" = "Trembling aspen",
                   "Black spruce" = "Black spruce",
                   "Minor species group" = "Minor species group")

newlabels2 <- list("changewithH" = expression(atop("Tree-level basal area growth rate",
                                                   paste("(", cm^{2}, " ", year^{-1}, " per tree)"))),
                   "mainTrend" = expression(atop("Plot-level basal area growth rate",
                                                 paste("(", cm^{2}, " ", year^{-1}, " per plot)"))))
figure_labeller <- function(variable,value){
  if(variable == "Species"){
    return(newlabels1[value])
  } else if (variable == "Direction"){
    return(newlabels2[value])
  }
}


figureA <- ggplot(data = output, aes(x = Year, y = PredictedBAGR))+
  facet_grid(Direction~Species, drop = FALSE, scales = "free_y", switch = "y",
             labeller = figure_labeller)+
  geom_ribbon(data = output[Direction == "mainTrend"],
              aes(x = Year, ymin = PredictedBAGR_Lower, ymax = PredictedBAGR_Upper),
              col = "white", fill = "gray", show.legend = FALSE)+
  geom_line(data = output[Direction == "mainTrend"],
            aes(x = Year, y = PredictedBAGR,
                linetype = as.factor(overallSignificant)),
            col = "black", size = 1, show.legend = FALSE)+
  geom_line(data = output[Direction != "mainTrend"], 
            aes(x = Year, y = PredictedBAGR, group = CompetitionIntensity, 
                col = CompetitionIntensity), 
            size = 1)+
  scale_color_continuous(name = "Competition \nintensity",
                         low = "#4d9221", high = "#c51b7d", 
                         breaks = c(3, 100),
                         labels = c("\nWeak", "Strong\n"),
                         guide = guide_colorbar(reverse = TRUE, ticks = FALSE))+
  scale_x_continuous(name = "Year", breaks = seq(1985, 2010, by = 5))+
  geom_segment(data = segmenttable, aes(x = x, xend = xend, y = y, yend = yend),
               size = 1)+
  geom_text(data = data.frame(Year = -Inf, PredictedBAGR = Inf,
                              Species = "All trees", texts = c("a", "b"),
                              Direction = c("changewithH", "mainTrend")),
            aes(x = Year, y = PredictedBAGR, label = texts), hjust = -1.5,
            vjust = 1, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_text(size = 15, face = "italic"),
        strip.text.y = element_text(size = 15),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.70, 0.85),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))




ggsave(file = file.path(workPath, "TablesFigures", 
                        paste("Figure S7. temporal trends in BA with H.png")),
       figureA,  width = 12, height = 10)







