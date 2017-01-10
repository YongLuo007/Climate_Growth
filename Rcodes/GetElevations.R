rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(raster)
workPath <- "~/GitHub/Climate_Growth"

masterTable <- read.csv(file.path(workPath, "data", "selectedPlotMasterTable.csv"),
                        header = TRUE, stringsAsFactors = FALSE) %>%
  data.table

source("~/GitHub/landwebNRV/landwebNRV/R/UTMtoLongLat.R")

locationData <- unique(masterTable[,.(PlotID, Easting, Northing, Zone = 14)], by = "PlotID")
locationData <- UTMtoLongLat(UTMTable = locationData)
locationData <- locationData$Transformed
demraster <- raster(file.path(workPath, "data", "StudyAreaClimates_BiomSIM", "MBDEM.tif"))
locationGeoPoint <- SpatialPointsDataFrame(coords = data.frame(locationData[,.(Longitude, Latitude)]),
                                           data = data.frame(locationData[,.(PlotID)]),
                                           proj4string = CRS("+proj=longlat"),
                                           match.ID = TRUE)

locationGeoPoint <- spTransform(locationGeoPoint, CRSobj = crs(demraster))
locationGeoPoint@data$elevations <- extract(demraster, locationGeoPoint)
locationData <- setkey(locationData, PlotID)[setkey(data.table(locationGeoPoint@data), PlotID),
                                                nomatch = 0]

locationData <- locationData[!is.na(elevations),.(Name = PlotID, ID = PlotID,
                                                  Latitude, Longitude,
                                                  Elevation = elevations)]

write.csv(locationData,
          file.path(workPath, "StudyAreaClimates_BiomSIM", "plotlocations.csv"),
          row.names = FALSE)




