rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(MuMIn)
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
source(file.path(workPath, "Rcodes", "Rfunctions","HeghyiCICalculation.R"))
inputData[PlotID == "46-132", Distance:=Distance/100]
CIcompetitionData <- inputData[!is.na(Distance),]
CIcompetitionData <- CIcompetitionData[,.(PlotNumber = PlotID, TreeNumber = TreeNumber,
                                          Year, Distance, Angle, DBH, Species, Biomass)]

  speciesData <- data.table::copy(analysesData)
  BiomassCIdata1_3 <- HeghyiCICalculation(data = CIcompetitionData,
                                       maxRadius = 12.62,
                                       sizeIndex = "Biomass",
                                       distanceWeight = 1.3)
  BiomassCIdata0_4 <- HeghyiCICalculation(data = CIcompetitionData,
                                          maxRadius = 12.62,
                                          sizeIndex = "Biomass",
                                          distanceWeight = 0.4)
  
  BiomassCIdata1_3[, uniTreeID:=paste(PlotNumber, "_", TreeNumber, sep = "")]
  BiomassCIdata0_4[, uniTreeID:=paste(PlotNumber, "_", TreeNumber, sep = "")]
  BiomassCIdata1_3 <- unique(BiomassCIdata1_3, by = c("Year", "uniTreeID"))[,InterH := NULL]
  BiomassCIdata0_4 <- unique(BiomassCIdata0_4, by = c("Year", "uniTreeID"))[,IntraH := NULL]
  BiomassCIdata <- setkey(BiomassCIdata0_4, Year, uniTreeID)[setkey(BiomassCIdata1_3, Year, uniTreeID), nomatch = 0]
  BiomassCIdata <- BiomassCIdata[,.(uniTreeID, IniYear = Year, IntraH, InterH, H = IntraH+InterH)]

  
  
  speciesData <- setkey(speciesData, uniTreeID, IniYear)[setkey(BiomassCIdata,
                                                                  uniTreeID, IniYear), 
                                                           nomatch = 0]

  

