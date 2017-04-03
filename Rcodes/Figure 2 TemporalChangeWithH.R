rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data",selectionMethod,
               "FullYearModels.RData"))
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
                     CompetitionIntensity = numeric(), PredictedABGR = numeric(),
                     PredictedABGR_Lower = numeric(), PredictedABGR_Upper = numeric(),
                     overallSignificant = numeric())
for(indispecies in studySpecies){
  speciesData <- analysesDataAll[Species == indispecies,]
  minABGR <- round(abs(min(speciesData$BiomassGR)), 3)+0.01
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
  mainTrends$PredictedABGR <- treeToPlotConvertor*(exp(fittedvalues$fit)-minABGR)
  mainTrends$PredictedABGR_Upper <- treeToPlotConvertor*(exp(fittedvalues$fit+1.98*fittedvalues$se.fit)-minABGR)
  mainTrends$PredictedABGR_Lower <- treeToPlotConvertor*(exp(fittedvalues$fit-1.98*fittedvalues$se.fit)-minABGR)
  if(speciesallHcoeff[rn == "Yearctd",]$`p-value` < 0.05){
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
    if(useLogQunatile95){
      H95quantilelower <- exp(as.numeric(quantile(log(speciesData$MidH), probs = 0.025)))
      H95quantileupper <- exp(as.numeric(quantile(log(speciesData$MidH), probs = 0.975)))
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
                                            H = exp(seq(log(min(speciesData$MidH)),
                                                        log(max(speciesData$MidH)),
                                                        length = 100)),
                                            stringsAsFactors = FALSE))
    }
    
    changewithH[,':='(Yearctd = Year-mean(speciesData$MidYear),
                      logDBHctd = 0,
                      logHctd = log(H)-mean(log(speciesData$MidH)),
                      logSActd = 0)]
    changewithH[,CompetitionIntensity:=as.numeric(as.factor(H))]
    fittedvalues <- predict(theallHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
    changewithH$PredictedABGR <- exp(fittedvalues$fit)-minABGR
    changewithH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)-minABGR
    changewithH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)-minABGR
    
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
output[,':='(Species = factor(Species, levels = studySpecies),
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

figure_labeller <- function(variable,value){
  if(variable == "Species"){
    return(newlabels1[value])
  } 
}


figureA <- ggplot(data = output[Direction == "mainTrend"], aes(x = Year, y = PredictedABGR))+
  facet_grid(.~Species, drop = TRUE,
             labeller = figure_labeller)+
  geom_ribbon(aes(x = Year, ymin = PredictedABGR_Lower, ymax = PredictedABGR_Upper),
              col = "white", fill = "gray", show.legend = FALSE)+
  geom_line(aes(x = Year, y = PredictedABGR,
                linetype = as.factor(overallSignificant)),
            col = "black", size = 1, show.legend = FALSE)+
   scale_y_continuous(name = expression(atop("Plot-level aboveground biomass growth rate",
                                             paste("(Kg ", year^{-1}, " per plot)"))),
                      limits = c(0, 82))+
  scale_x_continuous(name = "Year", breaks = seq(1985, 2010, by = 5))+
  geom_text(data = data.frame(Year = -Inf, PredictedABGR = Inf,
                              Species = "All trees", texts = "b"),
            aes(x = Year, y = PredictedABGR, label = texts), hjust = -1.5, vjust = 1, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 15),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_blank())




cutpoint <- 0.6
output_B <- output[Direction != "mainTrend"][
  , Panel:= factor("Upper", levels = c("Upper", "Lower"))]
output_BLower <- output_B[PredictedABGR < cutpoint,][, 
                        Panel:= factor("Lower", levels = c("Upper", "Lower"))]
output_B <- rbind(output_B, output_BLower)

newlabels1 <- list("All trees" = "All trees",
                   "Jack pine" = "Jack pine",
                   "Trembling aspen" = "Trembling aspen",
                   "Black spruce" = "Black spruce",
                   "Minor species group" = "Minor species group")
newlabels2 <- list("Upper" = "a",
                   "Lower" = "b")
figureB_labeller <- function(variable,value){
  if(variable == "Species"){
    return(newlabels1[value])
  } else if (variable == "Panel"){
    return(newlabels2[value])
  }
}

figureB <- ggplot(data = output_B[Panel == "Upper"], aes(x = Year, y = PredictedABGR))+
  facet_grid(Panel~Species, drop = TRUE,
             scale = "free_y",
             labeller = figureB_labeller)+
  geom_line(aes(group = CompetitionIntensity, col = CompetitionIntensity), 
            size = 1)+
  scale_y_continuous(name = expression(atop("Tree-level aboveground biomass growth rate",
                                            paste("(Kg ", year^{-1}, " per tree)"))))+
  scale_x_continuous(name = "Year", breaks = seq(1985, 2010, by = 5))+
  scale_color_continuous(name = "Competition \nintensity",
                         low = "#4d9221", high = "#c51b7d", 
                         breaks = c(3, 100),
                         labels = c("\nWeak", "Strong\n"),
                         guide = guide_colorbar(reverse = TRUE, ticks = FALSE))+
  geom_text(data = data.frame(x = -Inf, y = Inf, Species = "All trees",
                              Panel = "Upper", texts = "a"),
            aes(x = x, y = y, label = texts), vjust = 1, hjust = -1.5, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 15, face = "italic"),
        strip.background = element_blank(),
        # legend.direction = "horizontal",
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.70, 0.75),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

a_Grob <- ggplotGrob(figureA)
b_Grob <- ggplotGrob(figureB)

a_Grob$widths <- b_Grob$widths

dev(4)
clearPlot()
plotlayout <- rbind(c(1), c(1), c(2), c(2))
figure <- grid.arrange(b_Grob, a_Grob, 
                       layout_matrix = plotlayout)
workPath <- "~/GitHub/Climate_Growth"
if(useLogQunatile95){
  ggsave(file = file.path(workPath, "TablesFigures", 
                          paste("Figure 2. temporal trends with H_95Quantile", selectionMethod, ".png")),
         figure,  width = 12, height = 10)
} else {
  ggsave(file = file.path(workPath, "TablesFigures", 
                          paste("Figure 2. temporal trends with H", selectionMethod, ".png")),
         figure,  width = 11, height = 10)
  
}






