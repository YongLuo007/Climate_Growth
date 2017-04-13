
rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table); library(parallel)
library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}

selectionMethod <- "Year10Analyses"

analysesData <- fread(file.path(workPath, "data", selectionMethod, "finalData10.csv"))

analysesData <- analysesData[allCensusLiveTree == "yes",]

studySpecies <- c("All species", "Jack pine", "Trembling aspen", 
                  "Black spruce", "Minor species")
climates <- c("ATA", "GSTA", "NONGSTA",
              "ACMIA", "GSCMIA", 
              "ACO2A")
fullClimateModels <- list()

for(indispecies in studySpecies){
  speciesData <- analysesData[Species_Group == indispecies,]
  for(indiclimate in climates){
    speciesData$climate <- c(speciesData[, indiclimate, with = FALSE])
    minABGR <- round(abs(min(speciesData$BiomassGR)), 3)+0.01
    speciesData[,':='(logY = log(BiomassGR+minABGR), 
                      logDBHctd = log(MidDBH)-mean(log(MidDBH)), 
                      Climatectd = climate-mean(climate),
                      logHctd = log(H)-mean(log(H)),
                      logSActd = log(MidFA)-mean(log(MidFA)))]
    if(indispecies == "All species"){
      fullClimateModel <- lme(logY~logDBHctd+logSActd+Climatectd+logHctd+logDBHctd:logSActd+
                                logDBHctd:logHctd+logDBHctd:Climatectd+logSActd:logHctd+
                                logSActd:Climatectd+Climatectd:logHctd,
                              data = speciesData,
                              random = ~1|PlotID/uniTreeID, 
                              control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    } else {
      fullClimateModel <- lme(logY~logDBHctd+logSActd+Climatectd+logHctd+logDBHctd:logSActd+
                                logDBHctd:logHctd+logDBHctd:Climatectd+logSActd:logHctd+
                                logSActd:Climatectd+Climatectd:logHctd,
                              data = speciesData,
                              random = ~1|PlotID/uniTreeID, 
                              control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    }
    
    fullClimateModels[[paste(indispecies, "_", indiclimate, sep = "")]] <- fullClimateModel
    cat("Species", indispecies, "and climate", indiclimate, "is done. \n")
    rm(fullClimateModel)
  }
}
allHFixedCoeff <- lapply(fullClimateModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})

save.image(file.path(workPath, "data",
                     selectionMethod, 
                     "fullClimateModels.RData"))
