rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(MuMIn);library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- "J:/MBgrowth"
}
load(file.path(workPath, "data", "MBPSP.RData"))
names(allPSP) <- c("PlotID", "FMU", "TWP", "RGE", "PlotNumber",
                   "StandType", "StructureType", "VegeType",
                   "SStructure", "Moist", "YearEstabish", "Easting",
                   "Northing", "TreeNumber", "Species", "Distance",
                   "Angle", "DBH", "Height", "Status", "Class",
                   "TreeAge", "Year", "PlotSize")
inputData <- data.table::copy(allPSP)
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
inputData <- inputData[PlotID %in% unique(analysesData$PlotID),]

inputData[Species == "AS",species:="black spruce"]
# "BA"
inputData[Species == "BA",species:="balsam poplar"]
# "BF" 
inputData[Species == "BF",species:="balsam fir"]
# "BO"
inputData[Species == "BO",species:="white oak"] # bur oak to
# "BS"
inputData[Species == "BS",species:="black spruce"]
# "EC"
inputData[Species == "EC",species:="western redcedar"] # cedar
# "JP"
inputData[Species == "JP",species:="jack pine"]
# "MM"
inputData[Species == "MM",species:="silver maple"]
# "RP"
inputData[Species == "RP",species:="red pine"]
# "TA"
inputData[Species == "TA",species:="trembling aspen"]
# "TL"
inputData[Species == "TL",species:="tamarack larch"]
# "WB"
inputData[Species == "WB",species:="white birch"]
# "WE"
inputData[Species == "WE",species:="white elm"]
# "WS"
inputData[Species == "WS",species:="white spruce"]
# for unknownspecies
inputData[is.na(species), species:="unknown"]

inputData <- inputData[Status == 1 | Status == 0,]

source("~/GitHub/landwebNRV/landwebNRV/R/biomassCalculation.R")
inputData$Biomass <- biomassCalculation(species = inputData$species,
                                        DBH = inputData$DBH)
source(file.path(workPath, "Rcodes", "Rfunctions","HeghyiCICalculationModified.R"))
inputData[PlotID == "46-132", Distance:=Distance/100]
CIcompetitionData <- inputData[!is.na(Distance),]
CIcompetitionData <- CIcompetitionData[,.(PlotNumber = PlotID, TreeNumber = TreeNumber,
                                          Year, Distance, Angle, DBH, Species, Biomass)]
set(analysesData, , c("H_DBH", "IntraH_DBH", "InterH_DBH", 
                      "H_Biomass", "IntraH_Biomass", "InterH_Biomass", "IntraH1_3", 
                      "InterH0_4", "H"), NULL)

output <- data.table(Species = character(), sizeWeight = character(),
                     disweight = numeric(), IntraHAIC = numeric(), InterHAIC = numeric())
for(sizeWeight in seq(1, 2, by = 1)){
  for(disweight in seq(0, 2, by = 0.1)){
    processCIdata <- data.table::copy(CIcompetitionData)
    CIdata <- HeghyiCICalculation(data = processCIdata,
                                  maxRadius = 12.62,
                                  sizeIndex = "Biomass",
                                  distanceWeight = disweight, 
                                  sizeWeight = sizeWeight)
    analysesDataAll <- data.table::copy(analysesData)
    CIdata <- CIdata[, .(uniTreeID=paste(PlotNumber, "_", TreeNumber, sep = ""),
                         IniYear = Year,
                         IntraH, InterH)] %>%
      unique(., by = c("uniTreeID", "IniYear"))
    analysesDataAll <- setkey(analysesDataAll, uniTreeID, IniYear)[setkey(CIdata, uniTreeID, IniYear),
                                                                   nomatch = 0]
    
    for(indispecies in c("All", "JP", "TA", "BS", "Other")){
      if(indispecies == "All"){
        speciesData <- data.table::copy(analysesDataAll)
      } else if(indispecies == "Other"){
        speciesData <- analysesDataAll[!(Species %in% c("JP", "TA", "BS")),]
      } else {
        speciesData <- analysesDataAll[Species == indispecies,]
      }
      speciesData[,':='(logY = log(BiomassGR),
                        logIntraHctd = log(IntraH+1) - mean(log(IntraH+1)),
                        logInterHctd = log(InterH+1) - mean(log(InterH+1)))]
      IntraHM <- lme(logY~logIntraHctd,
                     random = ~1|PlotID,
                     data = speciesData,
                     control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
      InterHM <- lme(logY~logInterHctd,
                     random = ~1|PlotID,
                     data = speciesData,
                     control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
      aictable <- data.table(Species = indispecies, sizeWeight = sizeWeight, 
                             disweight = disweight, IntraHAIC = AIC(IntraHM), 
                             InterHAIC = AIC(InterHM))
      output <- rbind(output, aictable)
    }
    cat("disweight: ", disweight, "is done. \n")  
  }
}
  


# write.csv(output, file.path("~/GitHub/Climate_Growth/Results/bestWeightsbothsizeanddistance.csv"), row.names = F)

a <- melt(output, id.vars = c("Species", "sizeWeight", "disweight"), 
          measure.vars = c("IntraHAIC", "InterHAIC"),
          value.name = "Value")
a[,':='(variable = factor(variable, levels = c("IntraHAIC", "InterHAIC"),
                          labels = c("Intraspecific competition", "Interspecific competition")),
        Species = factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                         labels = c("All species", "Jack pine", "Trembling aspen", "Black spruce",
                                    "Other species")))]

a[,minvalue:=min(Value), by = c("Species", "variable")]

minvaluepoints <- a[Value == minvalue,]

figure <- ggplot(data=a, aes(x = disweight, y = Value))+
  facet_grid(Species~variable, scales = "free_y")+
  geom_line(aes(group = sizeWeight, col = sizeWeight), size = 1)+
  geom_point(data = minvaluepoints, aes(x = sizeWeight, y = Value, col = disweight),
             show.legend = FALSE, size = 2)+
  
  # geom_text(data = minvaluepoints[sizeIndex == "Biomass",], 
  #           aes(x = disweight, y = Value, label = disweight),
  #           vjust = -1, size = 5)+
  scale_y_continuous(name = "AIC")+
  scale_x_continuous(name = "Distance weight")+
  # scale_color_manual(name = "Size", values = c("black", "gray"))+
  guides(colour = guide_legend(title.position = "top", direction = "horizontal"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.position = c(0.1, 0.75),
        legend.background = element_rect(colour = "black"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"))

ggsave(file = file.path("~/GitHub/Climate_Growth/TablesFigures/bestdistanceweight.png"),
       figure, width = 11.5, height = 10)

