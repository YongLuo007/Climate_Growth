rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
if(!dir.exists(file.path(workPath, "data", selectionMethod))){
  dir.create(file.path(workPath, "data", selectionMethod))
}

firstRun <- FALSE
if(firstRun){
  analysesData <- fread(file.path(workPath, "data", "newAllDataRescaledComp1.csv"))
  analysesData[Species == "Other species", Species := "Minor species"]
  analysesData <- analysesData[allCensusLiveTree == "yes" & positiveGrowthTree == "yes",]
  myFunction <- function(sizeWeight, disWeight, analysesData, workPath){
    output <- data.table(Species = character(), sizeWeight = character(),
                         disWeight = numeric(), HAIC = numeric(), IntraHAIC = numeric(),
                         InterHAIC = numeric())
    for(i in sizeWeight){
      for(j in disWeight){
        newCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
                                     paste("CompetitionData_DW",j, "_SW", i, ".csv", sep = "")))
        analysesDataAll <- setkey(analysesData, uniTreeID, IniYear)[setkey(newCIdata, uniTreeID, IniYear),
                                                                    nomatch = 0]
        
        for(indispecies in c("All species", "Jack pine", "Trembling aspen",
                             "Black spruce", "Minor species")){
          speciesData <- analysesDataAll[Species == indispecies,]
          speciesData[,':='(logY = log(BiomassGR),
                            logHctd = log(H+1)-mean(log(H+1)),
                            logIntraHctd = log(IntraH+1) - mean(log(IntraH+1)),
                            logInterHctd = log(InterH+1) - mean(log(InterH+1)))]
          HModel <- lme(logY~logHctd,
                        random = ~1|uniTreeID,
                        data = speciesData,
                        control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
          IntraHModel <- lme(logY~logIntraHctd,
                             random = ~1|uniTreeID,
                             data = speciesData,
                             control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
          InterHModel <- lme(logY~logInterHctd,
                             random = ~1|uniTreeID,
                             data = speciesData,
                             control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
          aictable <- data.table(Species = indispecies, sizeWeight = i, 
                                 disWeight = j, HAIC = AIC(HModel),
                                 IntraHAIC = AIC(IntraHModel), 
                                 InterHAIC = AIC(InterHModel))
          output <- rbind(output, aictable)
        }
      }
    }
    return(output)
  }
  
  inputWeights <- list()
  m <- 1
  for(i in seq(0, 10, by = 0.1)){
    for(j in seq(0, 2, by = 0.1)){
      inputWeights[[m]] <- c(i, j)
      m <- m+1
    }
  }
  
  
  analysesData[,':='(H = NULL, IntraH = NULL, InterH = NULL)]
  
  
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c("lme", "AIC", "lmeControl", "data.table", "myFunction",
                                "setkey", "fread", "workPath", "analysesData"))
  
  allresults <- parLapply(cl, inputWeights, 
                          function(y) myFunction(sizeWeight = y[1], disWeight = y[2], 
                                                 analysesData = analysesData, workPath = workPath))
  stopCluster(cl)
  
  for(i in 1:length(allresults)){
    if(i == 1){
      output <- allresults[[i]]
    } else {
      output <- rbind(output, allresults[[i]])
    }
  }
  
  write.csv(output, file.path(workPath, "data", selectionMethod,
                              paste("bestAandB.csv", sep = "")),
            row.names = F)
} else {
  
  output <- fread(file.path(file.path(workPath, "data", selectionMethod,
                                      paste("bestAandB.csv", sep = ""))))
}


a <- melt(output, id.vars = c("Species", "sizeWeight", "disWeight"), 
          measure.vars = c("HAIC", "IntraHAIC", "InterHAIC"),
          value.name = "Value")

a[,minvalue:=min(Value), by = c("Species", "variable")]
bestWeightTable <- a[Value == minvalue,]

analysesData <- fread(file.path(workPath, "data", "newAllDataRescaledComp1.csv"))
studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce",
                  "Minor species")
analysesData[Species == "Other species", Species:="Minor species"]
analysesData[,':='(H = NULL, IntraH = NULL, InterH = NULL)]
for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies,]
  HsizeWeight <- bestWeightTable[Species == indispecies & variable == "HAIC",]$sizeWeight
  HdisWeight <- bestWeightTable[Species == indispecies & variable == "HAIC",]$disWeight
  HCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
                             paste("CompetitionData_DW",HdisWeight, "_SW", HsizeWeight, ".csv", sep = "")))
  
  #  IntraHsizeWeight <- bestWeightTable[Species == indispecies & variable == "IntraHAIC",]$sizeWeight
  #  IntraHdisWeight <- bestWeightTable[Species == indispecies & variable == "IntraHAIC",]$disWeight
  #  IntraHCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
  #                                  paste("CompetitionData_DW",IntraHdisWeight, "_SW", IntraHsizeWeight, ".csv", sep = "")))
  # # 
  #  InterHsizeWeight <- bestWeightTable[Species == indispecies & variable == "InterHAIC",]$sizeWeight
  #  InterHdisWeight <- bestWeightTable[Species == indispecies & variable == "InterHAIC",]$disWeight
  #  InterHCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
  #                                  paste("CompetitionData_DW",InterHdisWeight, "_SW", InterHsizeWeight, ".csv", sep = "")))
  #  allCIdata <- setkey(HCIdata[,.(uniTreeID, IniYear, H)], uniTreeID, IniYear)[setkey(IntraHCIdata[,.(uniTreeID, IniYear, IntraH)],
  #                                                                                    uniTreeID, IniYear), nomatch = 0]
  #  allCIdata <- setkey(allCIdata, uniTreeID, IniYear)[setkey(InterHCIdata[,.(uniTreeID, IniYear, InterH)],
  #                                                            uniTreeID, IniYear), nomatch = 0]
  # 
  speciesData <- setkey(speciesData, uniTreeID, IniYear)[setkey(HCIdata, uniTreeID, IniYear), 
                                                         nomatch = 0]
  if(indispecies == "All species"){
    newAllData <- speciesData
  } else {
    newAllData <- rbind(newAllData, speciesData)
  }
}

write.csv(newAllData, file.path(workPath, "data", selectionMethod, "finalData.csv"),
          row.names = FALSE)


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
  
  speciesData[,':='(logY = log(BiomassGR), 
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
                     "BestYearModels.RData"))


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
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
allHbestClimateModels <- list()
indiHbestClimateModels <- list()
allHbestIDVs <- list()
indiHbestIDVs <- list()

climates <- c("ATA", "GSTA", "NONGSTA",
              "ACMIA", "GSCMIA", 
              "ACO2A")
for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies,]
  for(indiclimate in climates){
    speciesData$climate <- c(speciesData[, indiclimate, with = FALSE])
    speciesData[,':='(logY = log(BiomassGR), 
                      logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                      Climatectd = climate-mean(climate),
                      logHctd = log(H+1)-mean(log(H+1)),
                      logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                      logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                      logSActd = log(IniFA+2.5)-mean(log(IniFA+2.5)))]
    allHoutput <- mixedModelSelection(DV = "logY", 
                                      IDV = c("logDBHctd", "Climatectd", "logHctd",
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
                                       IDV = c("logDBHctd", "Climatectd", "logIntraHctd",
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
    
    allHbestClimateModels[[paste(indispecies, "_", indiclimate, sep = "")]] <- allHbestModel
    indiHbestClimateModels[[paste(indispecies, "_", indiclimate, sep = "")]] <- indiHbestModel
    allHbestIDVs[[paste(indispecies, "_", indiclimate, sep = "")]] <- allHoutput$bestIDV
    indiHbestIDVs[[paste(indispecies, "_", indiclimate, sep = "")]] <- indiHoutput$bestIDV
    cat("Species", indispecies, "and climate", indiclimate, "is done. \n")
    rm(allHoutput, allHbestFormula, allHbestModel, indiHoutput, indiHbestFormu, indiHbestModel)
  }
}
allHFixedCoeff <- lapply(allHbestClimateModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})
indiHFixedCoeff <- lapply(indiHbestClimateModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})
save.image(file.path(workPath, "data", selectionMethod,
                     "BestClimateModels.RData"))



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


analysesData <- analysesData[allCensusLiveTree == "yes",]


studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Minor species")
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
for(indispecies in studySpecies){
  
  speciesData <- analysesData[Species == indispecies,]
  
  speciesData[,':='(logY = log(BiomassGR+15), 
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
  if(indispecies == "All species"){
    allHmodelSelectionSummaries <- allHoutput$modelSummary[, Species:=indispecies]
    allHbestFormulas <- list(allHbestFormula)
    names(allHbestFormulas) <- "All species"
    allHbestModels <- list(allHbestModel)
    names(allHbestModels) <- "All species"
  } else {
    allHmodelSelectionSummaries <- rbind(allHmodelSelectionSummaries, 
                                         allHoutput$modelSummary[, Species:=indispecies])
    allHbestFormulas[[indispecies]] <- allHbestFormula
    allHbestModels[[indispecies]] <- allHbestModel
  }
  cat("Species", indispecies, "is done. \n")
  rm(allHoutput, indiHoutput, allHbestFormula, allHbestModel, indiHbestFormu, indiHbestModel)
}

allHFixedCoeff <- lapply(allHbestModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})
rm(getIC, getICformFomula, indispecies, mixedModelSelection, 
   speciesData)
save.image(file.path(workPath, "data", selectionMethod,
                     "BestYearModels_IncludingNegativeGrowth.RData"))


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
climates <- c("ATA", "GSTA", "NONGSTA",
              "ACMIA", "GSCMIA", 
              "ACO2A")
allHbestClimateModels <- list()
indiHbestClimateModels <- list()
for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies,]
  for(indiclimate in climates){
    speciesData$climate <- c(speciesData[, indiclimate, with = FALSE])
    speciesData[,':='(logY = log(BiomassGR), 
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
    # allHbestIDVs[[paste(indispecies, "_", indiclimate, sep = "")]] <- allHoutput$bestIDV
    # indiHbestIDVs[[paste(indispecies, "_", indiclimate, sep = "")]] <- indiHoutput$bestIDV
    cat("Species", indispecies, "and climate", indiclimate, "is done. \n")
    rm(allHbestModel)
  }
}
allHFixedCoeff <- lapply(allHbestClimateModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})
indiHFixedCoeff <- lapply(indiHbestClimateModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})

save.image(file.path(workPath, "data",
                     "AllCensus_PositiveGrowth_RandomPlotADTree", 
                     "fullClimateModels.RData"))




