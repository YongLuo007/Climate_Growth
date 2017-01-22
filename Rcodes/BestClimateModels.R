rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table); library(parallel)
library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
analysesData <- read.csv(file.path(workPath, "data", "newAllDataRescaledComp1.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
analysesData <- analysesData[allCensusLiveTree == "yes" & positiveGrowthTree == "yes",]
studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Other species")
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
                                      random = ~1|uniTreeID, 
                                      control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    allHbestFormula <- as.formula(paste("logY~", paste(allHoutput$bestIDV, collapse = "+")))
    allHbestModel <- lme(fixed = allHbestFormula,
                         data = speciesData,
                         random = ~1|uniTreeID, 
                         control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    
    indiHoutput <- mixedModelSelection(DV = "logY", 
                                       IDV = c("logDBHctd", "Climatectd", "logIntraHctd",
                                               "logInterHctd",
                                               "logSActd"),
                                       maxInteraction = 2,
                                       ICTerm = "AIC",
                                       ICCut = 2,
                                       data = speciesData,
                                       random = ~1|uniTreeID, 
                                       control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    indiHbestFormu <- as.formula(paste("logY~", paste(indiHoutput$bestIDV, collapse = "+")))
    indiHbestModel <- lme(fixed = indiHbestFormu,
                          data = speciesData,
                          random = ~1|uniTreeID, 
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
save.image(file.path(workPath, "data", "ClimateModelSelection168Plots.RData"))
