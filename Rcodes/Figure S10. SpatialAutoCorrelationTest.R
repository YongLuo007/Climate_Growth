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
indispecies <- "All species"
allHbestFormula <- allHbestFormulas[[indispecies]]
theallHmodel <- allHbestModels[[indispecies]]
speciesData <- analysesData[Species == indispecies,]
speciesData[,':='(logY = log(BiomassGR), 
                  logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                  Yearctd = 0,
                  logHctd = log(H+1)-mean(log(H+1)),
                  logSActd = log(IniFA+2.5)-mean(log(IniFA+2.5)))]

speciesData$fittelogY <- predict(theallHmodel, newdata = speciesData,
                                 level = 0, se.fit = FALSE)
speciesData[, logYResiduals:=logY-fittelogY]
selectedplots <- unique(speciesData$PlotID)
for(indiplot in selectedplots){
  indiplotdata <- speciesData[PlotID == indiplot,]
  linearRegression <- lm(logYResiduals~Year, data = indiplotdata)
  plotcoeffs <- as.data.table(summary(linearRegression)$coefficients,
                              keep.rownames = TRUE)
  plotcoeffs <- plotcoeffs[rn == "Year",.(PlotID = indiplot, 
                                          slope = exp(Estimate)-1,
                                          P = `Pr(>|t|)`)]
  if(indiplot == selectedplots[1]){
    allslopes <- plotcoeffs
  } else {
    allslopes <- rbind(allslopes, plotcoeffs)
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
allslopeswithLocation <- setkey(locations, PlotID)[setkey(allslopes, PlotID),
                                                   nomatch = 0]
library(ade4)
plot_dists <- dist(cbind(allslopeswithLocation$Longitude,
                         allslopeswithLocation$Latitude))
slope_dists <- dist(allslopeswithLocation$slope)
a <- mantel.rtest(plot_dists, slope_dists, nrepet = 9999)
canadamap <- readRDS(file.path(workPath,
                               "data",
                               "canadamap.rds"))

MBlocation <- SpatialPoints(allslopeswithLocation[,.(Easting, Northing)],
                            proj4string = CRS("+proj=utm +zone=14 datum=NAD83"))
MBlocation <- spTransform(MBlocation, crs(canadamap))
MBlocation <- as.data.frame(MBlocation@coords)
MBlocation$slope <- allslopeswithLocation$slope

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

makeupG <- data.frame(expand.grid(y = seq(-4, -1, length = 20),
                                  slope = seq(min(allslopeswithLocation$slope),
                                              max(allslopeswithLocation$slope),
                                              length = 50)))

subfig <- ggplot(data=allslopeswithLocation)+
  geom_histogram(aes(slope, col = slope), bins = 20, 
                 col = "white", fill = "gray")+
  geom_raster(data=makeupG, aes(x=slope, y=y, fill=slope))+
  scale_fill_gradient2(low="red", mid = "gray", high = "blue", 
                       guide = FALSE)+
  geom_segment(data=data.frame(y = 0, yend = 0, slope = -Inf, slopeend = Inf),
               aes(x=slope, xend=slopeend, y=y, yend=yend))+
  geom_segment(data=data.frame(y = 0, yend = Inf, slope = -Inf, slopeend = -Inf),
               aes(x=slope, xend=slopeend, y=y, yend=yend))+
  scale_y_continuous(name = "Frequency", limits = c(-20, 35),
                     breaks = seq(0, 35, by = 7))+
  scale_x_continuous(name = expression(paste("Slope (kg ", year^{-2}, ")")),
                     breaks = seq(-0.1, 0.1, length = 5))+
  geom_text(data = data.frame(y = -8, slope = seq(-0.1, 0.1, length = 5),
                              texts = sprintf("%0.2f", round(seq(-0.1, 0.1, length = 5), 2))),
            aes(x = slope, y = y, label = texts), size = 5, angle = 45)+
  geom_point(data = data.frame(y = -4.5, slope = seq(-0.1, 0.1, length = 5)),
             aes(x = slope, y = y), pch = 17)+
  annotate("text", x = 0.02, y = -15, label = "Slope~(kg~year^{-2})", parse = TRUE,
           size = 5)+
  annotate("text", x = -0.11, y = -20, label = "Mantel~test:~italic(P)==0.97", parse = TRUE,
           size = 5, hjust = 0)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

subfigGrob <- ggplotGrob(subfig)
Figure1 <- ggplot(data = studyAreaall, aes(x = long, y = lat)) +
  geom_polygon(data = studyAreaall, aes(group = group), col = "gray", 
               alpha = 0, linetype = 1, size = 1)+
  geom_point(data = MBlocation,
             aes(x = Easting, y = Northing, col = slope), size = 2)+
  scale_color_gradient2(low="red", mid = "gray", high = "blue",
                        guide = FALSE)+
  geom_text(data = data.frame(y = Inf, x = -Inf, texts = "a"),
            aes(x = x, y = y, label = texts), size = 10, hjust = -1.5, vjust = 1.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

FigureA <- Figure1+
  annotation_custom(grob = subfigGrob, xmin = -20000, xmax = 410000,
                    ymin = 1600000, ymax = 2300000)

SAinfor <- analysesData[,.(minSA=min(IniFA), maxSA = max(FinFA)), by = PlotID]
SAinfor[,midSA:=(minSA+maxSA)/2]

allslopeswithLocation <- setkey(allslopeswithLocation, PlotID)[setkey(SAinfor, PlotID),
                                                               nomatch = 0]

lmcoef <- as.data.table(summary(lm(slope~minSA, data = allslopeswithLocation))$coefficients,
                        keep.rownames = TRUE)
lmfunction <- paste("slope==", round(lmcoef$Estimate[1], 4),"+", 
                    round(lmcoef$Estimate[2], 4), "~SA~(", "italic(P)==",
                    round(lmcoef$`Pr(>|t|)`[2], 2), ")", sep = "")
FigureB <- ggplot(data = allslopeswithLocation, aes(x = minSA, y = slope))+
  geom_point(col = "gray")+
  geom_smooth(method = "lm", col = "black")+
  annotate("text", x = 100, y = -0.05, label = lmfunction, parse = TRUE, size = 5)+
  geom_text(data = data.frame(y = Inf, x = -Inf, texts = "b"),
            aes(x = x, y = y, label = texts), size = 10, hjust = -1.5, vjust = 1.5)+
  scale_y_continuous(name = expression(paste("Slope (kg ", year^{-2}, ")")))+
  scale_x_continuous(name = "Stand age (SA, years)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

plotlayout <- cbind(1,2)
a <- grid.arrange(FigureA, FigureB, layout_matrix = plotlayout)
ggsave(file = file.path(workPath, "TablesFigures", 
                        paste("Figure S10. Spatial Autocorrelation check and SA check.png")),
       a, width = 12, height = 8)


