rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
inputData <- fread(file.path(workPath, "data", "forcompetitionIndex.csv"))
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
source(file.path(workPath, "Rcodes",  "Rfunctions", "biomassCalculation.R"))
inputData$Biomass <- biomassCalculation(species = inputData$species,
                                        DBH = inputData$DBH)
source(file.path(workPath, "Rcodes", "Rfunctions", "HeghyiCICalculationModified.R"))
inputData[PlotID == "46-132", Distance:=Distance/100]
CIcompetitionData <- inputData[!is.na(Distance),]
CIcompetitionData <- CIcompetitionData[,.(PlotNumber = PlotID, TreeNumber = TreeNumber,
                                          Year, Distance, Angle, DBH, Species, Biomass)]

sizeWeightList <- list("sizeWeight1" = seq(0, 2, by = 0.1),
                       "sizeWeight2" = seq(2.1, 4, by = 0.1),
                       "sizeWeight3" = seq(4.1, 6, by = 0.1),
                       "sizeWeight4" = seq(6.1, 8, by = 0.1),
                       "sizeWeight5" = seq(8.1, 10, by = 0.1))

for(k in 1:2){
  sizeWeight <-  sizeWeightList[[k]]
  disWeight  <-  seq(0, 2, by = 0.1)
  processCIdata <- data.table::copy(CIcompetitionData)
  CIdata <- HeghyiCICalculation(data = processCIdata,
                                maxRadius = 12.62,
                                sizeIndex = "Biomass",
                                distanceWeight = disWeight, 
                                sizeWeight = sizeWeight,
                                assymetricScale = "Rescale",
                                testing = FALSE)
  dd <- data.table::copy(CIdata$output)
  dd[, ':='(uniTreeID = paste(PlotNumber, "_", TreeNumber, sep = ""),
            IniYear = Year)][,':='(PlotNumber = NULL, TreeNumber = NULL, Year = NULL)]
  dd <- unique(dd, by = c("uniTreeID", "IniYear"))
  names(dd) <- gsub("\\.", "_", names(dd))
  
  for(i in sizeWeight){
    for(j in disWeight){
      newCIdata <- data.table::copy(dd)
      indicombweight <- paste(c("H_", "IntraH_", "InterH_"), paste("DW", j, "_SW", i, sep = ""), sep = "")
      indicombweight <- gsub("\\.", "_", indicombweight)
      setnames(newCIdata, indicombweight, c("H", "IntraH", "InterH"))
      newCIdata <- newCIdata[,.(uniTreeID, IniYear, H, IntraH, InterH)]
      write.csv(newCIdata, file.path(workPath, "data", "newAllCompetitioinData",
                                     paste("CompetitionData_DW", j, "_SW", i, ".csv", sep = "")))
    }
  }
}







