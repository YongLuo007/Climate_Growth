rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2, lib.loc = "~/GitHub/Climate_Growth/RequiredRPackages"); 
library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data", selectionMethod, "FullYearModels.RData"))
##### for overall temporal trends and its dependency on DBH and RBI
analysesDataAll[,SA:=IniFA+2.5]
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
  speciesData <- analysesDataAll[Species_Group == indispecies,]
  speciesallHcoeff <- allHcoeff[Species == indispecies, ]
  theallHmodel <- allfullModelsAll[[indispecies]]
  minBiomass <- round(abs(min(speciesData$BiomassGR)), 3)+0.01
  # change with H
  changewithH <- data.table(expand.grid(Species = indispecies,
                                        Direction = "changewithDBH",
                                        Year = seq(min(speciesData$IniYear), 
                                                   max(speciesData$FinYear),
                                                   length = 100), 
                                        DBH = exp(seq(log(min(speciesData$IniDBH)),
                                                      log(max(speciesData$FinDBH)),
                                                      length = 100)),
                                        stringsAsFactors = FALSE))
  changewithH[,':='(Yearctd = Year-mean(speciesData$MidYear),
                    logDBHctd = log(DBH)-mean(log(speciesData$MidDBH)),
                    logHctd = 0,
                    logSActd = 0)]
  changewithH[,CompetitionIntensity:=as.numeric(as.factor(DBH))]
  fittedvalues <- predict(theallHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
  changewithH$PredictedABGR <- exp(fittedvalues$fit)-minBiomass
  changewithH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)-minBiomass
  changewithH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)-minBiomass
  outputaddon <- changewithH[,.(Species, Direction, Year, 
                                CompetitionIntensity,
                                PredictedABGR, PredictedABGR_Lower,
                                PredictedABGR_Upper,
                                overallSignificant = 1)]
  if(speciesallHcoeff[Direction == "changewithDBH",]$`p-value`>=0.05){
    outputaddon[, overallSignificant := 2]
  }
  output <- rbind(output, outputaddon)
  
  rm(fittedvalues, changewithH, outputaddon)
  
  
  changewithH <- data.table(expand.grid(Species = indispecies,
                                        Direction = "changewithSA",
                                        Year = seq(min(speciesData$IniYear), 
                                                   max(speciesData$FinYear),
                                                   length = 100), 
                                        SA = exp(seq(log(min(speciesData$IniFA)),
                                                     log(max(speciesData$FinFA)),
                                                     length = 100)),
                                        stringsAsFactors = FALSE))
  changewithH[,':='(Yearctd = Year-mean(speciesData$MidYear),
                    logDBHctd = 0,
                    logHctd = 0,
                    logSActd = log(SA)-mean(log(speciesData$MidFA)))]
  changewithH[,CompetitionIntensity:=as.numeric(as.factor(SA))]
  fittedvalues <- predict(theallHmodel, newdata = changewithH, level = 0, se.fit = TRUE)
  changewithH$PredictedABGR <- exp(fittedvalues$fit)-minBiomass
  changewithH$PredictedABGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)-minBiomass
  changewithH$PredictedABGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)-minBiomass
  outputaddon <- changewithH[,.(Species, Direction, Year, 
                                CompetitionIntensity,
                                PredictedABGR, PredictedABGR_Lower,
                                PredictedABGR_Upper,
                                overallSignificant = 1)]
  if(speciesallHcoeff[Direction == "changewithSA",]$`p-value`>=0.05){
    outputaddon[, overallSignificant := 2]
  }
  output <- rbind(output, outputaddon)
  rm(fittedvalues, changewithH, outputaddon)
  
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

figure <- ggplot(data = output[overallSignificant == 1,], aes(x = Year, y = PredictedABGR))+
  facet_grid(Direction~Species, scale = "free_y",  drop = TRUE)+
  geom_line(aes(group = CompetitionIntensity, col = CompetitionIntensity), 
            size = 1)+
  scale_y_continuous(name = expression(atop("Tree-level aboveground biomass growth rate",
                                            paste("(Kg ", year^{-1}, " per tree)"))))+
  scale_x_continuous(name = "Year", breaks = seq(1985, 2010, by = 5))+
  scale_color_continuous(name = "DBH / SA",
                         low = "#4d9221", high = "#c51b7d",  breaks = c(3, 100),
                         labels = c("\nSmall / young", "Big / old\n"),
                         guide = guide_colorbar(reverse = TRUE, ticks = FALSE))+
  geom_segment(data = segmenttable, aes(x = x, xend = xend, y = y, yend = yend), size = 1)+
  geom_text(data = texttable, aes(x = x, y = y, label = texts), vjust = 1, size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 16),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_text(size = 15, face = "italic"),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.90, 0.25),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))
workPath <- "~/GitHub/Climate_Growth"
ggsave(file = file.path(workPath, "TablesFigures", 
                        paste("Figure S3. temporal trends with DBHSA", selectionMethod, ".png")),
       figure,  width = 11, height = 9.5)
