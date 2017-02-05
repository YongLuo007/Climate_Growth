rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra);library(grid)
library(raster);library(maptools);library(rgeos)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
load(file.path(workPath, "data",selectionMethod,
               "BestYearModels.RData"))
workPath <- "~/GitHub/Climate_Growth"
for(indispecies in studySpecies){
  allHbestFormula <- allHbestFormulas[[indispecies]]
  theallHmodel <- allHbestModels[[indispecies]]
  speciesData <- analysesData[Species == indispecies,]
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = 0,
                    logHctd = log(H)-mean(log(H)),
                    logSActd = log(IniFA+2.5)-mean(log(IniFA+2.5)))]
  speciesData$fittelogY <- predict(theallHmodel, newdata = speciesData,
                                   level = 0, se.fit = FALSE)
  plotRandom <- as.data.table(random.effects(theallHmodel)$PlotID,
                              keep.rownames = TRUE)
  names(plotRandom) <- c("PlotID", "PlotEffect")
  treeRandom <- as.data.table(random.effects(theallHmodel)$uniTreeID,
                              keep.rownames = TRUE)
  names(treeRandom) <- c("uniTreeID", "TreeEffect")
  treeRandom[, uniTreeID:=unlist(lapply(uniTreeID, 
                                  function(s) unlist(strsplit(s, split = "/", fixed = T))[2]))]
  speciesData <- setkey(speciesData, PlotID)[setkey(plotRandom, PlotID),
                                             nomatch = 0]
  speciesData <- setkey(speciesData, uniTreeID)[setkey(treeRandom, uniTreeID),
                                             nomatch = 0]
  
  speciesData[, logYResiduals:=logY-fittelogY-PlotEffect-TreeEffect]
  
  selectedplots <- unique(speciesData$PlotID)
  for(indiplot in selectedplots){
    indiplotdata <- speciesData[PlotID == indiplot,]
    linearRegression <- lm(logYResiduals~Year+logSActd:Year+Year:logHctd+Year:logDBHctd, data = indiplotdata)
    plotcoeffs <- as.data.table(summary(linearRegression)$coefficients,
                                keep.rownames = TRUE)
    plotcoeffs <- plotcoeffs[rn == "Year",.(Species = indispecies,
                                            PlotID = indiplot, 
                                            slope = Estimate,
                                            P = `Pr(>|t|)`)]
    if(indispecies == studySpecies[1] & indiplot == selectedplots[1]){
      allslopes <- plotcoeffs
    } else {
      allslopes <- rbind(allslopes, plotcoeffs)
    }
  }
}

locations <- read.csv(file.path(workPath, "data", "selectedPlotMasterTable.csv"),
                      header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table
locations <- locations[,.(PlotID, Zone = 14, Easting, Northing)] %>%
  unique(., by = "PlotID")
source(file.path(workPath, "Rcodes", "Rfunctions", "UTMtoLongLat.R"))
locations <- UTMtoLongLat(locations, "+proj=longlat")$Transformed
library(ade4)

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

allslopesOg <- data.table::copy(allslopes)
allslopes <- data.table::copy(allslopesOg)
makeupRaster <- data.frame(expand.grid(y = seq(1.5, 2.5, length = 50),
                                       x = seq(min(allslopes$slope),
                                               max(allslopes$slope),
                                               length = 100)))

legendPlot <- ggplot(data = data.frame(x = c(min(allslopes$slope), 
                                             max(allslopes$slope)),
                                       y = c(0, 5)), aes(x = x, y = y)) +
  geom_point(col = "white")+
  geom_rect(data = data.frame(x = -Inf, xmax = Inf, y = 0.7, ymax = 3.3), 
            aes(xmin = x, xmax = xmax, ymin = y, ymax = ymax), col = "black", fill = "white")+
  geom_raster(data = makeupRaster, aes(x = x, y = y, fill = x))+
  scale_fill_gradient2(low="red", mid = "gray", high = "green", guide = "none")+
  geom_point(data = data.frame(y = 1.3, x = seq(min(allslopes$slope), 
                                                max(allslopes$slope),
                                                length = 7),
                               texts = sprintf("%0.2f", 
                                               round(seq(min(allslopes$slope), 
                                                         max(allslopes$slope),
                                                         length = 7), 2))),
             aes(x = x, y = y), pch = 17, size = 2)+
  geom_text(data = data.frame(y = 1, x = seq(min(allslopes$slope), 
                                             max(allslopes$slope),
                                             length = 7),
                              texts = sprintf("%0.2f", 
                                              round(seq(min(allslopes$slope), 
                                                        max(allslopes$slope),
                                                        length = 7), 2))),
            aes(x = x, y = y, label = texts), size = 5)+
  geom_text(data = data.frame(y = 3, x = min(allslopes$slope),
                              texts = "Effect of Year on residuals"),
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
allfigures <- list()

for(indispecies in studySpecies){
  indispeciesslopes <- allslopesOg[Species==indispecies,][,Species:=NULL]
  allslopeswithLocation <- setkey(locations, PlotID)[setkey(indispeciesslopes, PlotID),
                                                     nomatch = 0]
  plot_dists <- dist(cbind(allslopeswithLocation$Longitude,
                           allslopeswithLocation$Latitude))
  slope_dists <- dist(allslopeswithLocation$slope)
  a <- mantel.rtest(plot_dists, slope_dists, nrepet = 9999)$pvalue
  mantelTestResults[[indispecies]] <- paste("Mantel~test:~italic(P)==", round(a, 2), sep = "")
  indispeciesLocations <- SpatialPoints(allslopeswithLocation[,.(Easting, Northing)],
                              proj4string = CRS("+proj=utm +zone=14 datum=NAD83"))
  indispeciesLocations <- spTransform(indispeciesLocations, crs(canadamap))
  indispeciesLocations <- as.data.frame(indispeciesLocations@coords) %>% data.table
  indispeciesLocations$slope <- allslopeswithLocation$slope
  indispeciesLocations$Species <- indispecies
  allfigures[[indispecies]] <- ggplot(data = studyAreaall, aes(x = long, y = lat)) +
    geom_polygon(data = studyAreaall, aes(group = group), col = "blue", 
                 alpha = 0, linetype = 2, size = 1)+
    geom_point(data = indispeciesLocations,
               aes(x = Easting, y = Northing, col = slope), size = 4)+
    scale_color_gradient2(low="red", mid = "gray", high = "green",
                          limits = c(min(allslopes$slope),
                                       max(allslopes$slope)),
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
a <- grid.arrange(allfigures[[1]], legendPlot, allfigures[[2]], allfigures[[3]],
                  allfigures[[4]],allfigures[[5]],
                  layout_matrix = plotlayout)
ggsave(file = file.path(workPath, "TablesFigures", 
                        paste("Figure S10. Spatial Autocorrelation check_Slope.png")),
       a, width = 12.5, height = 14.5)


