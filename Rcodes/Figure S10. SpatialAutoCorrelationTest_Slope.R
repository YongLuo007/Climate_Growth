rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra);library(grid)
library(raster);library(maptools);library(rgeos)
library(ade4)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data",selectionMethod,
               "FullYearModels.RData"))
workPath <- "~/GitHub/Climate_Growth"
locations <- read.csv(file.path(workPath, "data", "selectedPlotMasterTable.csv"),
                      header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table
locations <- locations[,.(PlotID, Zone = 14, Easting, Northing)] %>%
  unique(., by = "PlotID")
source(file.path(workPath, "Rcodes", "Rfunctions", "UTMtoLongLat.R"))
locations <- UTMtoLongLat(locations, "+proj=longlat")$Transformed


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



canadamap <- readRDS(file.path(workPath,
                               "data",
                               "canadamap.rds"))
provinces <- c("Manitoba")
ecozones <- readRDS(file.path(workPath, "data", "ecozones.rds"))
ecozoneses <- c("Boreal Cordillera", "Boreal PLain",  "Boreal Shield",
                "Taiga Shield", "Taiga Plain", "Taiga Cordillera",
                "Hudson Plain")
boreal <- ecozones[ecozones@data$ZONE_NAME %in% ecozoneses,]
boreal <- spTransform(boreal, CRSobj = crs(canadamap))
MBProvince <- canadamap[canadamap@data$NAME  == "Manitoba", ]
studyArea <- gIntersection(boreal, MBProvince)
studyAreaall <- fortify(studyArea, region = "id") %>% data.table


makeupRaster_PlotSlope <- data.frame(expand.grid(y = seq(1.5, 2.5, length = 50),
                                       x = seq(min(allRandom$PlotSlope),
                                               max(allRandom$PlotSlope),
                                               length = 100)))

legendPlot_PlotSlope <- ggplot(data = data.frame(x = c(min(allRandom$PlotSlope), 
                                             max(allRandom$PlotSlope)),
                                       y = c(0, 5)), aes(x = x, y = y)) +
  geom_point(col = "white")+
  geom_rect(data = data.frame(x = -Inf, xmax = Inf, y = 0.7, ymax = 3.3), 
            aes(xmin = x, xmax = xmax, ymin = y, ymax = ymax), col = "black", fill = "white")+
  geom_raster(data = makeupRaster_PlotSlope, aes(x = x, y = y, fill = x))+
  scale_fill_gradient2(low="red", mid = "gray", high = "green", guide = "none")+
  geom_point(data = data.frame(y = 1.3, x = seq(min(allRandom$PlotSlope), 
                                                max(allRandom$PlotSlope),
                                                length = 7),
                               texts = sprintf("%0.2f", 
                                               round(seq(min(allRandom$PlotSlope), 
                                                         max(allRandom$PlotSlope),
                                                         length = 7), 2))),
             aes(x = x, y = y), pch = 17, size = 2)+
  geom_text(data = data.frame(y = 1, x = seq(min(allRandom$PlotSlope), 
                                             max(allRandom$PlotSlope),
                                             length = 7),
                              texts = sprintf("%0.2f", 
                                              round(seq(min(allRandom$PlotSlope), 
                                                        max(allRandom$PlotSlope),
                                                        length = 7), 2))),
            aes(x = x, y = y, label = texts), size = 5)+
  geom_text(data = data.frame(y = 3, x = min(allRandom$PlotSlope),
                              texts = "Plot-level random slope of Year"),
            aes(x = x, y = y, label = texts), hjust = 0, size = 7)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())



k <- 1
mantelTestResults <- list()
allfigures_PlotSlope <- list()

for(indispecies in studySpecies){
  indispeciesslopes <- allRandom[Species==indispecies,][,Species:=NULL]
  allslopeswithLocation <- setkey(indispeciesslopes, PlotID)[setkey(locations, PlotID),
                                                     nomatch = 0]
  plot_dists <- dist(cbind(allslopeswithLocation$Longitude,
                           allslopeswithLocation$Latitude))
  slope_dists <- dist(allslopeswithLocation$PlotSlope)
  a <- mantel.rtest(plot_dists, slope_dists, nrepet = 9999)$pvalue
  mantelTestResults[[indispecies]] <- paste("Mantel~test:~italic(P)==", round(a, 2), sep = "")
  indispeciesLocations <- SpatialPoints(allslopeswithLocation[,.(Easting, Northing)],
                              proj4string = CRS("+proj=utm +zone=14 datum=NAD83"))
  indispeciesLocations <- spTransform(indispeciesLocations, crs(canadamap))
  indispeciesLocations <- as.data.frame(indispeciesLocations@coords) %>% data.table
  indispeciesLocations$slope <- allslopeswithLocation$PlotSlope
  indispeciesLocations$Species <- indispecies
  allfigures_PlotSlope[[indispecies]] <- ggplot(data = studyAreaall, aes(x = long, y = lat)) +
    geom_polygon(data = studyAreaall, aes(group = group), col = "blue", 
                 alpha = 0, linetype = 2, size = 1)+
    geom_point(data = indispeciesLocations,
               aes(x = Easting, y = Northing, col = slope), size = 4)+
    scale_color_gradient2(low="red", mid = "gray", high = "green",
                          limits = c(min(allRandom$PlotSlope),
                                       max(allRandom$PlotSlope)),
                          guide = "none")+
    geom_text(data = data.frame(y = Inf, x = -Inf, texts = letters[k]),
              aes(x = x, y = y, label = texts), size = 10, hjust = -1.5, vjust = 1.5)+
    annotate("text", x = Inf, y = -Inf, 
             label = mantelTestResults[[indispecies]], parse = TRUE,
             size = 5, vjust = -1.5, hjust = 1.2)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  k <- k+1
  if(indispecies == studySpecies[1]){
    allSpeciesLocations <- indispeciesLocations
  } else {
    allSpeciesLocations <- rbind(allSpeciesLocations, indispeciesLocations)
  }
}

plotlayout <- rbind(c(1,2), c(3, 4), c(5, 6))
allfiguresPlotSlope <- grid.arrange(allfigures_PlotSlope[[1]], legendPlot_PlotSlope, 
                                    allfigures_PlotSlope[[2]], allfigures_PlotSlope[[3]],
                  allfigures_PlotSlope[[4]],allfigures_PlotSlope[[5]],
                  layout_matrix = plotlayout)
ggsave(file = file.path(workPath, "TablesFigures", 
                        paste("Figure S10. Spatial Autocorrelation check_PlotSlope.png")),
       allfiguresPlotSlope, width = 12.5, height = 14.5)


makeupRaster_meanTreeSlope <- data.frame(expand.grid(y = seq(1.5, 2.5, length = 50),
                                                     x = seq(min(allRandom$MeanTreeSlope),
                                                             max(allRandom$MeanTreeSlope),
                                                             length = 100)))

legendPlot_meanTreeSlope <- ggplot(data = data.frame(x = c(min(allRandom$MeanTreeSlope), 
                                                           max(allRandom$MeanTreeSlope)),
                                                     y = c(0, 5)), aes(x = x, y = y)) +
  geom_point(col = "white")+
  geom_rect(data = data.frame(x = -Inf, xmax = Inf, y = 0.7, ymax = 3.3), 
            aes(xmin = x, xmax = xmax, ymin = y, ymax = ymax), col = "black", fill = "white")+
  geom_raster(data = makeupRaster_meanTreeSlope, aes(x = x, y = y, fill = x))+
  scale_fill_gradient2(low="red", mid = "gray", high = "green", guide = "none")+
  geom_point(data = data.frame(y = 1.3, x = seq(min(allRandom$MeanTreeSlope), 
                                                max(allRandom$MeanTreeSlope),
                                                length = 7),
                               texts = sprintf("%0.2f", 
                                               round(seq(min(allRandom$MeanTreeSlope), 
                                                         max(allRandom$MeanTreeSlope),
                                                         length = 7), 2))),
             aes(x = x, y = y), pch = 17, size = 2)+
  geom_text(data = data.frame(y = 1, x = seq(min(allRandom$MeanTreeSlope), 
                                             max(allRandom$MeanTreeSlope),
                                             length = 7),
                              texts = sprintf("%0.2f", 
                                              round(seq(min(allRandom$MeanTreeSlope), 
                                                        max(allRandom$MeanTreeSlope),
                                                        length = 7), 2))),
            aes(x = x, y = y, label = texts), size = 5)+
  geom_text(data = data.frame(y = 3, x = min(allRandom$MeanTreeSlope),
                              texts = "Average tree-level random slope of Year"),
            aes(x = x, y = y, label = texts), hjust = 0, size = 7)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())



k <- 1
mantelTestResults <- list()
allfigures_meanTreeSlope <- list()

for(indispecies in studySpecies){
  indispeciesslopes <- allRandom[Species==indispecies,][,Species:=NULL]
  allslopeswithLocation <- setkey(indispeciesslopes, PlotID)[setkey(locations, PlotID),
                                                             nomatch = 0]
  plot_dists <- dist(cbind(allslopeswithLocation$Longitude,
                           allslopeswithLocation$Latitude))
  slope_dists <- dist(allslopeswithLocation$MeanTreeSlope)
  a <- mantel.rtest(plot_dists, slope_dists, nrepet = 9999)$pvalue
  mantelTestResults[[indispecies]] <- paste("Mantel~test:~italic(P)==", round(a, 2), sep = "")
  indispeciesLocations <- SpatialPoints(allslopeswithLocation[,.(Easting, Northing)],
                                        proj4string = CRS("+proj=utm +zone=14 datum=NAD83"))
  indispeciesLocations <- spTransform(indispeciesLocations, crs(canadamap))
  indispeciesLocations <- as.data.frame(indispeciesLocations@coords) %>% data.table
  indispeciesLocations$slope <- allslopeswithLocation$MeanTreeSlope
  indispeciesLocations$Species <- indispecies
  allfigures_meanTreeSlope[[indispecies]] <- ggplot(data = studyAreaall, aes(x = long, y = lat)) +
    geom_polygon(data = studyAreaall, aes(group = group), col = "blue", 
                 alpha = 0, linetype = 2, size = 1)+
    geom_point(data = indispeciesLocations,
               aes(x = Easting, y = Northing, col = slope), size = 4)+
    scale_color_gradient2(low="red", mid = "gray", high = "green",
                          limits = c(min(allRandom$MeanTreeSlope),
                                     max(allRandom$MeanTreeSlope)),
                          guide = "none")+
    geom_text(data = data.frame(y = Inf, x = -Inf, texts = letters[k]),
              aes(x = x, y = y, label = texts), size = 10, hjust = -1.5, vjust = 1.5)+
    annotate("text", x = Inf, y = -Inf, 
             label = mantelTestResults[[indispecies]], parse = TRUE,
             size = 5, vjust = -1.5, hjust = 1.2)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  k <- k+1
  if(indispecies == studySpecies[1]){
    allSpeciesLocations <- indispeciesLocations
  } else {
    allSpeciesLocations <- rbind(allSpeciesLocations, indispeciesLocations)
  }
}

plotlayout <- rbind(c(1,2), c(3, 4), c(5, 6))
allfiguresMeanTreeSlope <- grid.arrange(allfigures_meanTreeSlope[[1]], legendPlot_meanTreeSlope, 
                                        allfigures_meanTreeSlope[[2]], allfigures_meanTreeSlope[[3]],
                                        allfigures_meanTreeSlope[[4]],allfigures_meanTreeSlope[[5]],
                                        layout_matrix = plotlayout)
ggsave(file = file.path(workPath, "TablesFigures", 
                        paste("Figure S11. Spatial Autocorrelation check_MeanTreeSlope.png")),
       allfiguresMeanTreeSlope, width = 12.5, height = 14.5)



