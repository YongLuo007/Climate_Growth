rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(raster)
library(maptools); library(rgeos); library(colorRamps)

workPath <- "~/GitHub/Climate_Growth"


canadamap <- readRDS(file.path(workPath,
                               "data",
                               "canadamap.rds"))
provinces <- c("Quebec", "Nova Scotia", "Saskatchewan", "Alberta", "Newfoundland and Labrador", 
               "British Columbia", "New Brunswick", "Prince Edward Island", 
               "Manitoba", "Ontario")
canadamap <- canadamap[canadamap@data$NAME %in% provinces,]
ecozones <- readRDS(file.path(workPath, "data", "ecozones.rds"))
ecozoneses <- c("Boreal Cordillera", "Boreal PLain",  "Boreal Shield",
                "Taiga Shield", "Taiga Plain", "Taiga Cordillera",
                "Hudson Plain")
boreal <- ecozones[ecozones@data$ZONE_NAME %in% ecozoneses,]
boreal <- spTransform(boreal, CRSobj = crs(canadamap))
boreal <- gIntersection(boreal, canadamap)



df<- data.frame(id = "boreal")
row.names(df) <- 1
boreal <- SpatialPolygonsDataFrame(boreal, data = df)


locations <- read.csv(file.path(workPath, "data", "selectedPlotMasterTable.csv"),
                      header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table

Year10Plots <- fread(file.path(workPath, "data", "Year10Analyses", "finalData10.csv"))
locations <- locations[PlotID %in% unique(Year10Plots$PlotID),]


locations <- locations[,.(PlotID, Easting, Northing)] %>%
  unique(., by = "PlotID")
MBlocation <- SpatialPoints(locations[,.(Easting, Northing)], proj4string = CRS("+proj=utm +zone=14 datum=NAD83"))
MBlocation <- spTransform(MBlocation, crs(canadamap))
studyAreaSquare <- crop(canadamap, MBlocation)

MBlocation <- as.data.frame(MBlocation@coords)

studyArea <- gIntersection(boreal, studyAreaSquare)

lines <- list()
i <- 1
for(longid in seq(-140, -20, by = 10)){
  lines[[i]] <- Lines(list(Line(coords = data.frame(x = longid, y = seq(20, 75, by = 0.05)))), 
                      ID = paste("long", longid, sep = ""))
  names(lines)[i] <- paste("long", longid, sep = "")
  i <- i+1
}
for(latid in seq(30, 70, 5)){
  lines[[i]] <- Lines(list(Line(coords = data.frame(x = seq(-140, -20, by = 0.5), y = latid))), 
                      ID = paste("lat", latid, sep = ""))
  names(lines)[i] <- paste("lat", latid, sep = "")
  i <- i+1
}
coordReference <- SpatialLines(lines, proj4string = CRS("+proj=longlat"))
coordReference <- spTransform(coordReference, crs(canadamap))
rm(i)
for(i in 1:length(coordReference)){
  a <- as.data.frame(coordReference@lines[[i]]@Lines[[1]]@coords)
  a$line <- names(lines)[i]
  if(i == 1){
    coordReferenceFinal <- a
  } else {
    coordReferenceFinal <- rbind(coordReferenceFinal, a)
  }
}
coordReferenceFinal <- data.table(coordReferenceFinal)
names(coordReferenceFinal)[1:2] <- c("long", "lat")

canadamapall <- fortify(canadamap, region = "NAME") %>%
  data.table
canadamapall[,fill:=1]
borealall <- fortify(boreal, region = "id") %>%
  data.table
studyAreaall <- fortify(studyArea, region = "id") %>% data.table

coordReferenceFinal <- data.table(coordReferenceFinal)[long <= max(canadamapall$long) &
                                                       long >= min(canadamapall$long) &
                                                       lat <= max(canadamapall$lat) &
                                                       lat >= min(canadamapall$lat),]

coordReferenceFinal <- coordReferenceFinal[!(line %in% c("long-140", "long-50", "lat40", "lat65")),]


# studyareaall <- fortify(studyAreaPoly, region = "CMIRate") 
coordReferenceLabels <- coordReferenceFinal[!(line %in% paste("long", seq(-90, -60,  by = 10), sep = "")),
                                            .(long = max(long), lat = max(lat)), by = line]
coordReferenceLabels <- rbind(coordReferenceLabels,
                              coordReferenceFinal[line %in% paste("long", seq(-90, -60, by = 10), sep = ""),][
                                ,.(long = min(long), lat = max(lat)), by = line])
coordReferenceLabels[line %in% paste("long", seq(-130, -60, by = 10), sep = ""),
                     labels := paste(seq(-130, -60, by = 10), "W", sep = "")]
coordReferenceLabels[line %in% paste("lat", seq(45, 60, by = 5), sep = ""),
                     labels := paste(seq(45, 60, by = 5), "N", sep = "")]

Figure1a <- ggplot(data = canadamapall, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group, fill = as.factor(fill)))+
  scale_fill_manual(values = "gray")+
  geom_polygon(data = borealall, aes(x = long, y = lat, group = group),
               fill = "deepskyblue")+
  geom_path(aes(group = group), col = "white", size = 1)+
  geom_line(data = coordReferenceFinal, aes(x = long, y = lat, group = line),
            colour = "gray", linetype = 2)+
  # geom_point(data = MBlocation,
  #            aes(x = MBlocation$Easting, y = MBlocation$Northing), col = "blue",
  #            pch = 1, size = 1.5)+
  geom_polygon(data = studyAreaall, aes(x = long, y = lat, group = group), col = "red", 
               alpha = 0, linetype = 2, size = 1)+
  geom_text(data = coordReferenceLabels, aes(x = long, y = lat, label = labels))+
  annotate("text", label = "BC", x = -1800000, y = 2000000, size = 5)+
  annotate("text", label = "AB", x = -1200000, y = 1700000, size = 5)+
  annotate("text", label = "SK", x = -680000, y = 1500000, size = 5)+
  annotate("text", label = "MB", x = -80000, y = 1500000, size = 5)+
  annotate("text", label = "ON", x = 600000, y = 1200000, size = 5)+
  annotate("text", label = "QC", x = 1500000, y = 1700000, size = 5)+
  annotate("text", label = "USA", x = -680000, y = 900000, size = 5)+
  annotate("text", label = "a", x = -Inf, y = Inf, size = 8, vjust = 1.5, hjust = -1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

monthlyClimate <- read.csv(file.path(workPath, "data", "background_monthlyclimates.csv"), header = TRUE,
                           stringsAsFactors = FALSE) %>% data.table
names(monthlyClimate)[5:9] <- c("YearMonth", "Temperature", "Prep", "snowfall", "PET")
monthlyClimate[,':='(Year = as.numeric(unlist(lapply(YearMonth, function(s) {unlist(strsplit(s, "/", fixed = T))[1]}))),
                     Month = as.numeric(unlist(lapply(YearMonth, function(s) {unlist(strsplit(s, "/", fixed = T))[2]}))))]
monthlyClimate[Month>=10, Year:=Year+1]
annualClimate <- monthlyClimate[!(Year %in% c(1980, 2013)),.(Latitude = mean(Latitude), Longitude = mean(Longitude),
                                   AnnualT = mean(Temperature), AnnualCMI = sum(Prep-PET)),
                                by = c("ID", "Year")]
IDs <- unique(annualClimate$ID)
output <- data.table(ID = numeric(), Latitude = numeric(), Longitude = numeric(),
                     TemperatureTrend = numeric(), CMITrend = numeric())
for(id in IDs){
  output <- rbind(output, 
                  data.table(ID = id, Latitude = annualClimate[ID == id,]$Latitude[1],
                             Longitude = annualClimate[ID == id,]$Longitude[1],
                             TemperatureTrend = coef(lm(AnnualT~Year, data = annualClimate[ID == id,]))[2],
                             CMITrend = coef(lm(AnnualCMI~Year, data = annualClimate[ID == id,]))[2]))
}
AnnualTRaster <- rasterFromXYZ(xyz = output[,.(Longitude, Latitude, TemperatureTrend)], crs = crs("+proj=longlat"))
AnnualTRaster <- projectRaster(AnnualTRaster,
                               res = c(0.05, 0.05),
                               crs = crs(AnnualTRaster))
AnnualTRaster <- projectRaster(AnnualTRaster, crs = crs(canadamap))
AnnualTRaster <- mask(AnnualTRaster, studyArea)

AnnualCMIRaster <- rasterFromXYZ(xyz = output[,.(Longitude, Latitude, CMITrend)], crs = crs("+proj=longlat"))
AnnualCMIRaster <- projectRaster(AnnualCMIRaster,
                               res = c(0.05, 0.05),
                               crs = crs(AnnualCMIRaster))
AnnualCMIRaster <- projectRaster(AnnualCMIRaster, crs = crs(canadamap))
AnnualCMIRaster <- mask(AnnualCMIRaster, studyArea)

AnnualTRaster_p <- rasterToPoints(AnnualTRaster)
AnnualTRaster_p <- data.frame(AnnualTRaster_p)
colnames(AnnualTRaster_p) <- c("long", "lat", "Temperature")
cutpoints1 <- c(min(AnnualTRaster_p$Temperature), seq(0.005, 0.065, by = 0.01), max(AnnualTRaster_p$Temperature))
AnnualTRaster_p$TemperatureClass <- cut(AnnualTRaster_p$Temperature, breaks = cutpoints1,
                                        labels = seq(0, 0.07, by = 0.01), include.lowest = TRUE)


AnnualCMIRaster_p <- rasterToPoints(AnnualCMIRaster)
AnnualCMIRaster_p <- data.frame(AnnualCMIRaster_p)
colnames(AnnualCMIRaster_p) <- c("long", "lat", "CMI")
cutpoints2 <- seq(min(AnnualCMIRaster_p$CMI), max(AnnualCMIRaster_p$CMI), length = 9)
labels2 <- c()
for(i in 1:8){
  labels2 <- c(labels2, (cutpoints2[i]+cutpoints2[i+1])/2)
}
AnnualCMIRaster_p$CMIClass <- cut(AnnualCMIRaster_p$CMI, breaks = cutpoints2,
                                  labels = round(labels2, 2), include.lowest = TRUE)

Figure1b <- ggplot(data = AnnualTRaster_p, aes(x = long, y = lat))+
  geom_raster(aes(fill = Temperature))+
  scale_fill_gradient2(name = expression(atop("Annual temperature trend ", 
                                         paste("(", degree, "C yea",r^{-1}, ")"))),
                       low = "#2166ac", high = "#b2182b")+
  guides(fill = guide_colourbar(title.position = "top", barwidth = unit(2, "in")))+
  geom_point(data = MBlocation,
             aes(x = MBlocation$Easting, y = MBlocation$Northing), col = "black",
             pch = 3, size = 2)+
  annotate("text", label = "b", x = -Inf, y = Inf, size = 8, vjust = 1.5, hjust = -1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(angle = 45, size = 12),
        legend.direction = "horizontal",
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.30, 0.13))

Figure1c <- ggplot(data = AnnualCMIRaster_p, aes(x = long, y = lat))+
  geom_raster(aes(fill = CMI))+
  scale_fill_gradient2(name = expression(atop("Annual climate moisture index trend  ", 
                                              paste("(mm yea",r^{-1}, ")"))),
                       low = "#b2182b", high = "#2166ac")+
  guides(fill = guide_colourbar(title.position = "top", barwidth = unit(2, "in")))+
  geom_point(data = MBlocation,
             aes(x = MBlocation$Easting, y = MBlocation$Northing), col = "black",
             pch = 3, size = 2)+
  annotate("text", label = "c", x = -Inf, y = Inf, size = 8, vjust = 1.5, hjust = -1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.30, 0.13))

library(gridExtra)
plotlayout <- rbind(c(1, 1), c(1, 1), c(2, 3), c(2, 3), c(2, 3))
dev(4)
a <- grid.arrange(Figure1a, Figure1b, Figure1c, layout_matrix = plotlayout)



ggsave(file = file.path(workPath, "TablesFigures", "Figure 1.png"), 
       a, width = 10, height = 10)









