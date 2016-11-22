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
                    logIntraHctd = log(IntraH1_3+1)-mean(log(IntraH1_3+1)),
                    logInterHctd = log(InterH0_4+1)-mean(log(InterH0_4+1)),
                    RBIctd = RBI - mean(RBI),
                    logSActd = log(SA)-mean(log(SA)),
                    logSBctd = log(PlotBiomass)-mean(log(PlotBiomass)))]
  if(indispecies == "All"){
    themodel <- lme(logY~logDBHctd+Yearctd+logIntraHctd+logInterHctd+RBIctd+
                      logDBHctd:logInterHctd+logDBHctd:RBIctd+
                      Yearctd:RBIctd+
                      logIntraHctd:RBIctd+logIntraHctd:RBIctd+
                      Yearctd:RBIctd:logIntraHctd+
                      Yearctd:RBIctd:logInterHctd+
                      logDBHctd:logIntraHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=20000, msMaxIter = 20000),
                    data = speciesData)
  } else if(indispecies == "JP"){
    themodel <- lme(logY~logDBHctd+Yearctd+logIntraHctd+logInterHctd+RBIctd+
                      logDBHctd:logIntraHctd+logDBHctd:RBIctd+
                      Yearctd:RBIctd+logIntraHctd:RBIctd+
                      Yearctd:RBIctd:logIntraHctd+
                      logDBHctd:logIntraHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000),
                    data = speciesData)
  } else if(indispecies == "TA"){
    themodel <- lme(logY~logDBHctd+Yearctd+logIntraHctd+logInterHctd+RBIctd+
                      logDBHctd:logInterHctd+logDBHctd:RBIctd+
                      Yearctd:RBIctd+logIntraHctd:RBIctd+
                      Yearctd:RBIctd:logIntraHctd+
                      logDBHctd:logIntraHctd:RBIctd+
                      logDBHctd:logInterHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000),
                    data = speciesData)
  } else if(indispecies == "BS"){
    themodel <- lme(logY~logDBHctd+logIntraHctd+logInterHctd+RBIctd+
                      logDBHctd:Yearctd+logDBHctd:RBIctd+
                      Yearctd:logIntraHctd+Yearctd:RBIctd+
                      logIntraHctd:RBIctd+logInterHctd:RBIctd+
                      Yearctd:RBIctd:logIntraHctd+
                      logDBHctd:logIntraHctd:RBIctd+
                      logDBHctd:logInterHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000),
                    data = speciesData)
  } else {
    themodel <- lme(logY~logDBHctd+Yearctd+logIntraHctd+logInterHctd+RBIctd+
                      logDBHctd:RBIctd+logInterHctd:RBIctd,
                    random = ~1+Yearctd|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000),
                    data = speciesData)
  }
  allbestmodels[[indispecies]] <- themodel
}

save.image(file.path(workPath, "Results", "bestYearModel.RData"))

