rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
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
source(file.path(workPath, "Rcodes",  "biomassCalculation.R"))
inputData$Biomass <- biomassCalculation(species = inputData$species,
                                        DBH = inputData$DBH)
source(file.path(workPath, "Rcodes", "Rfunctions", "HeghyiCICalculationModified.R"))
inputData[PlotID == "46-132", Distance:=Distance/100]
CIcompetitionData <- inputData[!is.na(Distance),]
CIcompetitionData <- CIcompetitionData[,.(PlotNumber = PlotID, TreeNumber = TreeNumber,
                                          Year, Distance, Angle, DBH, Species, Biomass)]
set(analysesData, , c("H_DBH", "IntraH_DBH", "InterH_DBH", 
                      "H_Biomass", "IntraH_Biomass", "InterH_Biomass", "IntraH1_3", 
                      "InterH0_4", "H"), NULL)
firstRun <- TRUE
if(firstRun){
  sizeWeight <-  seq(5.1, 8, by = 0.1)
  disweight  <-  seq(0, 2, by = 0.1)
  processCIdata <- data.table::copy(CIcompetitionData)
  CIdata <- HeghyiCICalculation(data = processCIdata,
                                maxRadius = 12.62,
                                sizeIndex = "Biomass",
                                distanceWeight = disweight, 
                                sizeWeight = sizeWeight,
                                assymetricScale = "Rescale")
  dd <- data.table::copy(CIdata)
  
  dd[, ':='(uniTreeID = paste(PlotNumber, "_", TreeNumber, sep = ""),
            IniYear = Year)][,':='(PlotNumber = NULL, TreeNumber = NULL, Year = NULL)]
  dd <- unique(dd, by = c("uniTreeID", "IniYear"))
  names(dd) <- gsub("\\.", "_", names(dd))
  
  output <- data.table(Species = character(), sizeWeight = character(),
                       disweight = numeric(), IntraHAIC = numeric(), InterHAIC = numeric())
  for(i in sizeWeight){
    for(j in disweight){
      newCIdata <- data.table::copy(dd)
      analysesDataAll <- data.table::copy(analysesData)
      indicombweight <- paste(c("H_", "IntraH_", "InterH_"), paste("DW", j, "_SW", i, sep = ""), sep = "")
      indicombweight <- gsub("\\.", "_", indicombweight)
      
      setnames(newCIdata, indicombweight, c("H", "IntraH", "InterH"))
      newCIdata <- newCIdata[,.(uniTreeID, IniYear, H, IntraH, InterH)]
      analysesDataAll <- setkey(analysesDataAll, uniTreeID, IniYear)[setkey(newCIdata, uniTreeID, IniYear),
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
        aictable <- data.table(Species = indispecies, sizeWeight = i, 
                               disweight = j, IntraHAIC = AIC(IntraHM), 
                               InterHAIC = AIC(InterHM))
        output <- rbind(output, aictable)
      }
      cat("sizeWeight: ", i, "disweight: ", j, "is done. \n")  
    }
  }
  
  write.csv(output, file.path(workPath, "bestWeightsbothsizeanddistanceRescaled2.csv"), row.names = F)
  rm(i, j, output)
  
}
 outputadd <- read.csv(file.path(workPath, "bestWeightsbothsizeanddistanceRescaled.csv"),
                   header = TRUE, stringsAsFactors = FALSE) %>% data.table
 output <- read.csv(file.path(workPath, "bestWeightsbothsizeanddistanceRescaled2.csv"),
                    header = TRUE, stringsAsFactors = FALSE) %>% data.table
 output <- rbind(outputadd, output) 
 write.csv(output, file.path(workPath, "bestWeightsbothsizeanddistanceRescaledAll.csv"), row.names = F)
 rm(i, j)
 
output <- unique(output, by = c("Species", "sizeWeight", "disweight"))
a <- melt(output, id.vars = c("Species", "sizeWeight", "disweight"), 
          measure.vars = c("IntraHAIC", "InterHAIC"),
          value.name = "Value")
a[,':='(variable = factor(variable, levels = c("IntraHAIC", "InterHAIC"),
                          labels = c("Intraspecific competition", "Interspecific competition")),
        Species = factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                         labels = c("All species", "Jack pine", "Trembling aspen", "Black spruce",
                                    "Other species")))]

a[,minvalue:=min(Value), by = c("Species", "variable")]

bestWeightTable <- a[Value == minvalue,]

for(indispecies in c("All species", "Jack pine", "Trembling aspen", "Black spruce",
                     "Other species")){
  processCIdata <- data.table::copy(CIcompetitionData)
  analysesDataAll <- data.table::copy(analysesData)
  
  bestIntraSizeWei <- bestWeightTable[Species == indispecies & 
                                        variable == "Intraspecific competition", ]$sizeWeight
  bestIntraDisWei <- bestWeightTable[Species == indispecies & 
                                       variable == "Intraspecific competition", ]$disweight
  bestInterSizeWei <- bestWeightTable[Species == indispecies & 
                                        variable == "Interspecific competition", ]$sizeWeight
  
  bestInterDisWei <- bestWeightTable[Species == indispecies & 
                                       variable == "Interspecific competition", ]$disweight
  bestIntraHdata <- HeghyiCICalculation(data = processCIdata,
                                        maxRadius = 12.62,
                                        sizeIndex = "Biomass",
                                        distanceWeight = bestIntraDisWei, 
                                        sizeWeight = as.numeric(bestIntraSizeWei),
                                        assymetricScale = "Rescale")
  bestInterHdata <- HeghyiCICalculation(data = processCIdata,
                                        maxRadius = 12.62,
                                        sizeIndex = "Biomass",
                                        distanceWeight = bestInterDisWei, 
                                        sizeWeight = as.numeric(bestInterSizeWei),
                                        assymetricScale = "Rescale")
  
  bestIntraHdata[, ':='(uniTreeID=paste(PlotNumber, "_", TreeNumber, sep = ""),
                        IniYear = Year)]
  bestIntraHdata[,':='(PlotNumber = NULL, TreeNumber = NULL, Year = NULL)]
  names(bestIntraHdata)[1:3] <- c("H", "IntraH", "InterH")
  bestIntraHdata <- unique(bestIntraHdata, by = c("uniTreeID", "IniYear"))
  bestInterHdata[, ':='(uniTreeID=paste(PlotNumber, "_", TreeNumber, sep = ""),
                        IniYear = Year)]
  bestInterHdata[,':='(PlotNumber = NULL, TreeNumber = NULL, Year = NULL)]
  names(bestInterHdata)[1:3] <- c("H", "IntraH", "InterH")
  bestInterHdata <- unique(bestInterHdata, by = c("uniTreeID", "IniYear"))
  bestIntraHdata[,':='(InterH = NULL, H = NULL)]
  bestInterHdata[,':='(IntraH = NULL, H = NULL)]
  allcidata <- setkey(bestIntraHdata, uniTreeID, IniYear)[setkey(bestInterHdata, uniTreeID, IniYear),
                                                          nomatch = 0]
  if(indispecies == "All species"){
    speciesData <- data.table::copy(analysesDataAll)
  } else if(indispecies == "Other species"){
    speciesData <- analysesDataAll[!(Species %in% c("JP", "TA", "BS")),]
  } else if (indispecies == "Jack pine"){
    speciesData <- analysesDataAll[Species == "JP",]
  } else if(indispecies == "Trembling aspen"){
    speciesData <- analysesDataAll[Species == "TA",]
  } else {
    speciesData <- analysesDataAll[Species == "BS",]
  }
  speciesData <- setkey(speciesData, uniTreeID, IniYear)[setkey(allcidata, uniTreeID, IniYear),nomatch = 0]
  speciesData <- cbind(data.table(DataType = rep(indispecies, nrow(speciesData))), speciesData)
  if(indispecies == "All species"){
    allnewdata <- speciesData
  } else {
    allnewdata <- rbind(allnewdata, speciesData)
  }
}
write.csv(allnewdata, 
          file.path(workPath, "newAllDataRescaledComp.csv"),
          row.names = F)
