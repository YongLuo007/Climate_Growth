rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table);library(MuMIn)
library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"

analysesData <- fread(file.path(workPath, "data", selectionMethod, "finalData.csv"))


analysesData <- analysesData[allCensusLiveTree == "yes" & positiveGrowthTree == "yes",]


studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Minor species")
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
for(indispecies in studySpecies){
  
  speciesData <- analysesData[Species == indispecies,]
  
  speciesData[,':='(logY = log(BAGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = Year-mean(Year),
                    logHctd = log(H+1)-mean(log(H+1)),
                    logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                    logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                    logSActd = log(IniFA+2.5)-mean(log(IniFA+2.5)))]
  allHoutput <- mixedModelSelection(DV = "logY", 
                                    IDV = c("logDBHctd", "Yearctd", "logHctd",
                                            "logSActd"),
                                    maxInteraction = 2,
                                    ICTerm = "AIC",
                                    ICCut = 2,
                                    data = speciesData,
                                    random = ~1|PlotID/uniTreeID, 
                                    control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  allHbestFormula <- as.formula(paste("logY~", paste(allHoutput$bestIDV, collapse = "+")))
  allHbestModel <- lme(fixed = allHbestFormula,
                       data = speciesData,
                       random = ~1|PlotID/uniTreeID, 
                       control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  
  indiHoutput <- mixedModelSelection(DV = "logY", 
                                     IDV = c("logDBHctd", "Yearctd", "logIntraHctd",
                                             "logInterHctd",
                                             "logSActd"),
                                     maxInteraction = 2,
                                     ICTerm = "AIC",
                                     ICCut = 2,
                                     data = speciesData,
                                     random = ~1|PlotID/uniTreeID, 
                                     control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  indiHbestFormu <- as.formula(paste("logY~", paste(indiHoutput$bestIDV, collapse = "+")))
  indiHbestModel <- lme(fixed = indiHbestFormu,
                        data = speciesData,
                        random = ~1|PlotID/uniTreeID, 
                        control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  if(indispecies == "All species"){
    allHmodelSelectionSummaries <- allHoutput$modelSummary[, Species:=indispecies]
    allHbestFormulas <- list(allHbestFormula)
    names(allHbestFormulas) <- "All species"
    allHbestModels <- list(allHbestModel)
    names(allHbestModels) <- "All species"
    
    indiHmodelSelectionSummaries <- indiHoutput$modelSummary[, Species:=indispecies]
    indiHbestFormulas <- list(indiHbestFormu)
    names(indiHbestFormulas) <- "All species"
    indiHbestModels <- list(indiHbestModel)
    names(indiHbestModels) <- "All species"
    
  } else {
    allHmodelSelectionSummaries <- rbind(allHmodelSelectionSummaries, 
                                         allHoutput$modelSummary[, Species:=indispecies])
    allHbestFormulas[[indispecies]] <- allHbestFormula
    allHbestModels[[indispecies]] <- allHbestModel
    
    indiHmodelSelectionSummaries <- rbind(indiHmodelSelectionSummaries, 
                                          indiHoutput$modelSummary[, Species:=indispecies])
    indiHbestFormulas[[indispecies]] <- indiHbestFormu
    indiHbestModels[[indispecies]] <- indiHbestModel
  }
  cat("Species", indispecies, "is done. \n")
  rm(allHoutput, indiHoutput, allHbestFormula, allHbestModel, indiHbestFormu, indiHbestModel)
}

allHFixedCoeff <- lapply(allHbestModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})
indiHFixedCoeff <- lapply(indiHbestModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})

rm(getIC, getICformFomula, indispecies, mixedModelSelection, 
   speciesData)
save.image(file.path(workPath, "data", selectionMethod,
                     "BestYearModels_BA.RData"))


rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table); library(parallel)
library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
analysesData <- fread(file.path(workPath, "data", selectionMethod, "finalData.csv"))
analysesData <- analysesData[allCensusLiveTree == "yes" & positiveGrowthTree == "yes",]
studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Minor species")
allHbestClimateModels <- list()
indiHbestClimateModels <- list()
climates <- c("ATA", "GSTA", "NONGSTA",
              "ACMIA", "GSCMIA", 
              "ACO2A")
for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies,]
  for(indiclimate in climates){
    speciesData$climate <- c(speciesData[, indiclimate, with = FALSE])
    speciesData[,':='(logY = log(BAGR), 
                      logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                      Climatectd = climate-mean(climate),
                      logHctd = log(H+1)-mean(log(H+1)),
                      logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                      logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                      logSActd = log(IniFA+2.5)-mean(log(IniFA+2.5)))]
    allHbestModel <- lme(logY~logDBHctd+logSActd+Climatectd+logHctd+logDBHctd:logSActd+
                           logDBHctd:logHctd+logDBHctd:Climatectd+logSActd:logHctd+
                           logSActd:Climatectd+Climatectd:logHctd,
                         
                         data = speciesData,
                         random = ~1|PlotID/uniTreeID, 
                         control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    
    indiHbestModel <- lme(logY~logDBHctd+logSActd+Climatectd+logIntraHctd+logInterHctd+
                            logDBHctd:logSActd+logDBHctd:Climatectd+logDBHctd:logIntraHctd+
                            logDBHctd:logInterHctd+logSActd:Climatectd+logSActd:logIntraHctd+
                            logSActd:logInterHctd+Climatectd:logIntraHctd+Climatectd:logInterHctd+
                            logIntraHctd:logInterHctd,
                          data = speciesData,
                          random = ~1|PlotID/uniTreeID, 
                          control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    
    allHbestClimateModels[[paste(indispecies, "_", indiclimate, sep = "")]] <- allHbestModel
    indiHbestClimateModels[[paste(indispecies, "_", indiclimate, sep = "")]] <- indiHbestModel
    cat("Species", indispecies, "and climate", indiclimate, "is done. \n")
    rm(allHbestModel, indiHbestModel)
  }
}
allHFixedCoeff <- lapply(allHbestClimateModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})
indiHFixedCoeff <- lapply(indiHbestClimateModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})
save.image(file.path(workPath, "data", selectionMethod,
                     "fullClimateModels_BA.RData"))





