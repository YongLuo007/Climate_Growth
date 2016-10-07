rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
studySpecies <- c("All", "JP", "BS", "TA", "Other")
allbestmodels <- list()

for(indispecies in studySpecies){
  if(indispecies == "All"){
    speciesData <- data.table::copy(analysesData)
  } else if(indispecies == "Other"){
    speciesData <- analysesData[!(Species %in% c("JP", "TA", "BS")),]
  } else {
    speciesData <- analysesData[Species == indispecies,]
  }
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = Year-mean(Year),
                    logHctd = log(Hegyi)-mean(log(Hegyi)),
                    RBIctd = RBI - mean(RBI))]
  if(indispecies == "All"){
    themodel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+
                      logDBHctd:Yearctd+logDBHctd:logHctd+
                      logDBHctd:RBIctd+Yearctd:logHctd+
                      Yearctd:RBIctd+logHctd:RBIctd+
                      logDBHctd:Yearctd:logHctd+
                      logDBHctd:Yearctd:RBIctd+
                      logDBHctd:logHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=20000, msMaxIter = 20000),
                    data = speciesData)
  } else if(indispecies == "JP"){
    themodel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+
                      logDBHctd:logHctd+logDBHctd:RBIctd+
                      Yearctd:RBIctd+logHctd:RBIctd+Yearctd:logHctd+
                      logDBHctd:Yearctd:logHctd+
                      logDBHctd:Yearctd:RBIctd+
                      logDBHctd:logHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000),
                    data = speciesData)
  } else if(indispecies == "TA"){
    themodel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+
                      logDBHctd:logHctd+logDBHctd:RBIctd+
                      Yearctd:RBIctd+logHctd:RBIctd+Yearctd:logHctd+
                      logDBHctd:logHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000),
                    data = speciesData)
  } else if(indispecies == "BS"){
    themodel <- lme(logY~logDBHctd+logHctd+RBIctd+logDBHctd:Yearctd+
                      logDBHctd:RBIctd+logHctd:Yearctd+RBIctd:Yearctd+
                      logHctd:RBIctd+logDBHctd:logHctd:Yearctd+
                      logDBHctd:RBIctd:Yearctd+logDBHctd:logHctd:RBIctd+
                      logHctd:RBIctd:Yearctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000),
                    data = speciesData)
  } else {
    themodel <- lme(logY~logDBHctd+Yearctd+logHctd+RBIctd+logDBHctd:RBIctd+
                      logHctd:RBIctd+logDBHctd:logHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000),
                    data = speciesData)
  }
  allbestmodels[[indispecies]] <- themodel
}

save.image(file.path(workPath, "Results", "bestYearModel.RData"))

