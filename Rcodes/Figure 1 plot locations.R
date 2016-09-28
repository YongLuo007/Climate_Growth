rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(raster)
library(maptools); library(rgeos)

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


locations <- read.csv(file.path(workPath, "data", "masterTable.csv"),
                      header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table
workingData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"),
                        header = TRUE,
                        stringsAsFactors = FALSE) %>%
  data.table
locations <- locations[PlotID %in% unique(workingData$PlotID),]

locations <- locations[,.(PlotID, Easting, Northing)] %>%
  unique(., by = "PlotID")
MBlocation <- SpatialPoints(locations[,.(Easting, Northing)], proj4string = CRS("+proj=utm +zone=14 datum=NAD83"))
MBlocation <- spTransform(MBlocation, crs(canadamap))
MBlocation <- as.data.frame(MBlocation@coords)

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

figure1 <- ggplot(data = canadamapall, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group, fill = as.factor(fill)))+
  scale_fill_manual(values = "gray")+
  geom_polygon(data = borealall, aes(x = long, y = lat, group = group),
               fill = "green")+
  geom_path(aes(group = group), col = "white", size = 1)+
  geom_line(data = coordReferenceFinal, aes(x = long, y = lat, group = line),
            colour = "gray", linetype = 2)+
  geom_point(data = MBlocation,
             aes(x = MBlocation$Easting, y = MBlocation$Northing), col = "blue",
             pch = 1, size = 2)+
  geom_text(data = coordReferenceLabels, aes(x = long, y = lat, label = labels))+
  annotate("text", label = "BC", x = -1800000, y = 2000000, size = 5)+
  annotate("text", label = "AB", x = -1200000, y = 1700000, size = 5)+
  annotate("text", label = "SK", x = -680000, y = 1500000, size = 5)+
  annotate("text", label = "MB", x = -80000, y = 1500000, size = 5)+
  annotate("text", label = "ON", x = 600000, y = 1200000, size = 5)+
  annotate("text", label = "QC", x = 1500000, y = 1700000, size = 5)+
  annotate("text", label = "USA", x = -680000, y = 900000, size = 5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

ggsave(file = file.path(workPath, "TablesFigures", "Figure 1.png"), figure1, width = 14, height = 6)









