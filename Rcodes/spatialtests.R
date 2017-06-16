rm(list = ls())
library(nlme)
library(data.table)
library(dplyr)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data",selectionMethod,
               "FullYearModels.RData"))


allresiduals <- lapply(allfullModelsAll, function(s) summary(s)$data[, 
                                                     Residuals := as.numeric(residuals(s))])
plotslope <- list()
treeslope <- list()
for(indispecies in studySpecies){
  if(indispecies == "All species"){
    indimodel <- allfullModelsAll[[indispecies]]
    indiplotslope <- data.table(random.effects(indimodel)$PlotID,
                                keep.rownames = TRUE)
  
    indiplotslope[, PlotID := unlist(lapply(rn, function(s) unlist(strsplit(s, "/", fixed = TRUE))[2]))]
    inditreeslope <- data.table(random.effects(indimodel)$uniTreeID,
                                keep.rownames = TRUE)
    inditreeslope[, PlotID := unlist(lapply(rn, function(s) unlist(strsplit(s, "/", fixed = TRUE))[2]))]
    plotslope[[indispecies]] <- indiplotslope
    treeslope[[indispecies]] <- inditreeslope
  } else {
    
    indimodel <- allfullModelsAll[[indispecies]]
    indiplotslope <- data.table(random.effects(indimodel)$PlotID,
                                keep.rownames = TRUE)
    
    indiplotslope[, PlotID := unlist(lapply(rn, function(s) unlist(strsplit(s, "/", fixed = TRUE))[1]))]
    inditreeslope <- data.table(random.effects(indimodel)$uniTreeID,
                                keep.rownames = TRUE)
    inditreeslope[, PlotID := unlist(lapply(rn, function(s) unlist(strsplit(s, "/", fixed = TRUE))[1]))]
    plotslope[[indispecies]] <- indiplotslope
    treeslope[[indispecies]] <- inditreeslope
  }
}


