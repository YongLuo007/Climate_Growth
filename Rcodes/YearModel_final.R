rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
load(file.path(workPath, "data", "theBestModels.RData"))
rm(allbestmodels)
allmodels <- list()

studySpecies <- c("All", "JP", "BS", "TA", "Other")

for(indispecies in studySpecies){
  if(indispecies == "All"){
    speciesData <- analysesData
  } else if(indispecies == "Other"){
    speciesData <- analysesData[!(Species %in% c("JP", "BS", "TA")),]
  } else {
    speciesData <- analysesData[Species == indispecies,]
  }
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = Year-mean(Year),
                    logHctd = log(Hegyi)-mean(log(Hegyi)),
                    RBIctd = RBI - mean(RBI))]
  
  fullmodel <- FALSE
  if(fullmodel){
    theModel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+
                      logDBHctd:Yearctd+logDBHctd:logHctd+
                      logDBHctd:RBIctd+Yearctd:logHctd+
                      Yearctd:RBIctd+logHctd:RBIctd+
                      logDBHctd:Yearctd:logHctd+
                      logDBHctd:logHctd:RBIctd+
                      Yearctd:logHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID,
                    data = speciesData,
                    control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  } else {
    if(indispecies == "All"){
      theModel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+
                        logDBHctd:logHctd+
                        logDBHctd:RBIctd+Yearctd:logHctd+
                        Yearctd:RBIctd+logHctd:RBIctd+
                        logDBHctd:Yearctd:logHctd+
                        logDBHctd:logHctd:RBIctd+
                        Yearctd:logHctd:RBIctd,
                      random = ~1+Yearctd|PlotID/uniTreeID,
                      data = speciesData,
                      control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
      
    } else if(indispecies == "JP"){
      theModel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+
                        logDBHctd:logHctd+
                        logDBHctd:RBIctd+
                        Yearctd:RBIctd+logHctd:RBIctd+
                        logDBHctd:Yearctd:logHctd+
                        logDBHctd:logHctd:RBIctd+
                        Yearctd:logHctd:RBIctd,
                      random = ~1+Yearctd|PlotID/uniTreeID,
                      data = speciesData,
                      control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    } else if(indispecies == "TA"){
      theModel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+
                        logDBHctd:logHctd+
                        logDBHctd:RBIctd+
                        Yearctd:RBIctd+logHctd:RBIctd+
                        logDBHctd:Yearctd:logHctd+
                        logDBHctd:logHctd:RBIctd+
                        Yearctd:logHctd:RBIctd,
                      random = ~1+Yearctd|PlotID/uniTreeID,
                      data = speciesData,
                      control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
      
    } else if(indispecies == "BS"){
      theModel <- lme(logY~logDBHctd+logHctd+RBIctd+
                        logDBHctd:RBIctd+Yearctd:logHctd+
                        Yearctd:RBIctd+logHctd:RBIctd+
                        logDBHctd:Yearctd:logHctd+
                        logDBHctd:logHctd:RBIctd,
                      random = ~1+Yearctd|PlotID/uniTreeID,
                      data = speciesData,
                      control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    } else if(indispecies == "Other"){
      theModel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+
                        logDBHctd:RBIctd+
                        Yearctd:RBIctd+logHctd:RBIctd+
                        logDBHctd:logHctd:RBIctd,
                      random = ~1+Yearctd|PlotID/uniTreeID,
                      data = speciesData,
                      control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    }
  }
  
  allmodels[[indispecies]] <- theModel
  
}

save.image(file.path(workPath, "data", "BiomassGR_Year_Model_Final.RData"))