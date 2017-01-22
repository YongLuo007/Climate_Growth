rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table);library(MuMIn)
library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
analysesData <- read.csv(file.path(workPath, "data", "finalData.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table

analysesData[,TreeMT:=length(IniYear), by = c("Species", "uniTreeID")]
analysesData <- analysesData[TreeMT>=2, ]
analysesData <- analysesData[positiveGrowthTree == "yes",]

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
                                    random = ~1|uniTreeID, 
                                    control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  allHbestFormula <- as.formula(paste("logY~", paste(allHoutput$bestIDV, collapse = "+")))
  allHbestModel <- lme(fixed = allHbestFormula,
                       data = speciesData,
                       random = ~1|uniTreeID, 
                       control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  
  indiHoutput <- mixedModelSelection(DV = "logY", 
                                     IDV = c("logDBHctd", "Yearctd", "logIntraHctd",
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
save.image(file.path(workPath, "data", "finalYearModels_RandomUniTreeID_TwoCensusPositive.RData"))


