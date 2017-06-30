rm(list = ls())
library(nlme)
library(data.table)
library(dplyr)
library(sp)
library(ade4)

workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data",selectionMethod,
               "FullYearModels.RData"))


allresiduals <- lapply(allfullModelsAll, function(s) summary(s)$data[, 
                                                     Residuals := as.numeric(residuals(s))])


locations <- read.csv(file.path(workPath, "data", "selectedPlotMasterTable.csv"),
                      header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table
locations <- locations[,.(PlotID, Zone = 14, Easting, Northing)] %>%
  unique(., by = "PlotID")
source(file.path(workPath, "Rcodes", "Rfunctions", "UTMtoLongLat.R"))
locations <- UTMtoLongLat(locations, "+proj=longlat")$Transformed[,.(PlotID, Longitude, Latitude)]

mantelTestResults <- data.table(Species = character(), P = numeric())

for(indispecies in studySpecies){
  speciesdata <- allresiduals[[indispecies]]
  speciesdata <- speciesdata[,.(meanRedisuals = mean(Residuals)), by = PlotID]
  speciesdata <- setkey(speciesdata, PlotID)[setkey(locations, PlotID), nomatch = 0]
  plot_dists <- dist(cbind(speciesdata$Longitude,
                           speciesdata$Latitude))
  Redisual_dists <- dist(speciesdata$meanRedisuals)
  a <- suppressWarnings(mantel.rtest(plot_dists, Redisual_dists, nrepet = 9999)$pvalue)
  mantelTestResults <- rbind(mantelTestResults,
                             data.table(Species = indispecies,
                                        P = round(as.numeric(a), 2)))
}


