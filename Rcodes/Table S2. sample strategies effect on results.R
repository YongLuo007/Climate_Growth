rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra);library(grid)
library(raster);library(maptools);library(rgeos)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data", selectionMethod, "FullYearModels.RData"))
workPath <- "~/GitHub/Climate_Growth"
MeasureInfor <- analysesDataAll[,.(MeasurementLength = max(FinYear)-min(IniYear),
                                     minYear = min(IniYear)),
                                  by = c("Species_Group", "PlotID")]
setnames(MeasureInfor, "Species_Group", "Species")
for(indispecies in studySpecies){
  theallHmodel <- allfullModelsAll[[indispecies]]
  plotRandom <- as.data.table(random.effects(theallHmodel)$PlotID,
                              keep.rownames = TRUE)
  names(plotRandom) <- c("PlotID", "PlotIntercept", "PlotSlope")
  plotRandom <- plotRandom[,.(Species = indispecies,
                              PlotID, PlotIntercept, PlotSlope)]
  treeRandom <- as.data.table(random.effects(theallHmodel)$uniTreeID,
                              keep.rownames = TRUE)
  names(treeRandom) <- c("PlotIDuniTreeID", "TreeIntercept", "TreeSlope")
  treeRandom[, PlotID:=unlist(lapply(PlotIDuniTreeID, 
                                     function(s) unlist(strsplit(s, split = "/", fixed = T))[1]))]
  treeRandom[, uniTreeID:=unlist(lapply(PlotIDuniTreeID, 
                                        function(s) unlist(strsplit(s, split = "/", fixed = T))[2]))]
  MeanTreeRandom <- treeRandom[,.(Species = indispecies,
                                  MeanTreeIntercept = mean(TreeIntercept),
                                  MeanTreeSlope = mean(TreeSlope)),
                               by = "PlotID"]
  indispeciesrandom <- setkey(plotRandom, PlotID, Species)[setkey(MeanTreeRandom, PlotID, Species),
                                                           nomatch = 0]
  if(indispecies == "All species"){
    allRandom <- indispeciesrandom
  } else {
    allRandom <- rbind(allRandom, indispeciesrandom)
  }
}
allslopes <- setkey(allRandom, Species, PlotID)[setkey(MeasureInfor, Species, PlotID),
                                                nomatch = 0]





for(indispecies in studySpecies){
  minYearModel_Year <- as.data.table(summary(lm(PlotSlope~minYear, 
                                                data = allslopes[Species == indispecies]))$coefficients,
                                     keep.rownames = TRUE)
  MeasurementLModel_Year <- as.data.table(summary(lm(PlotSlope~MeasurementLength, 
                                              data = allslopes[Species == indispecies]))$coefficients,
                                   keep.rownames = TRUE)
  
  YearEffect <- rbind(minYearModel_Year, MeasurementLModel_Year)
  YearEffect <- YearEffect[rn %in% c("NofTree", "minYear", "minSA", "MeasurementLength"),][,.(Variable = rn,
                                                                         YearCorrelation = paste(round(Estimate, 3),
                                                                                             "(P = ",round(`Pr(>|t|)`, 4),")",
                                                                                             sep = ""))]
  minYearModel_Interaction <- as.data.table(summary(lm(MeanTreeSlope~minYear, 
                                                data = allslopes[Species == indispecies]))$coefficients,
                                     keep.rownames = TRUE)
  MeasurementLModel_Interaction <- as.data.table(summary(lm(MeanTreeSlope~MeasurementLength, 
                                                     data = allslopes[Species == indispecies]))$coefficients,
                                          keep.rownames = TRUE)
  
  Interaction <- rbind(minYearModel_Interaction,
                       MeasurementLModel_Interaction)
  Interaction <- Interaction[rn %in% c("NofTree", "minYear", "minSA", "MeasurementLength"),][,.(Variable = rn,
                                                                           InteractionCorrelation = paste(round(Estimate*(10^3), 2),
                                                                                                          "*10^-3 (P = ",round(`Pr(>|t|)`, 4),")",
                                                                                                          sep = ""))]
  alleffects <- setkey(YearEffect, Variable)[setkey(Interaction, Variable),
                                             nomatch = 0]
  alleffects[,Species:=indispecies]
  if(indispecies == "All species"){
    alloutput <- alleffects
  } else {
    alloutput <- rbind(alloutput, alleffects)
  }

}

write.csv(alloutput, file.path(workPath, "TablesFigures",
                               "SampleStrategyEffects.csv"),
          row.names = FALSE)


