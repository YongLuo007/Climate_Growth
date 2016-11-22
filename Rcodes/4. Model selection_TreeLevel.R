rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table);library(MuMIn)
library(parallel)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "newAllData.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Other species")
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
for(indispecies in studySpecies){
  speciesData <- analysesData[DataType == indispecies,]
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    RBIctd = RBI - mean(RBI),
                    Yearctd = Year-mean(Year),
                    logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                    logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                    logSActd = log(SA)-mean(log(SA)),
                    logSBctd = log(PlotBiomass)-mean(log(PlotBiomass)))]
  tempoutput <- mixedModelSelection(DV = "logY", 
                                    IDV = c("logDBHctd", "Yearctd", "logIntraHctd",
                                            "logInterHctd", "logSActd", "RBIctd"),
                                    maxInteraction = 2,
                                    ICTerm = "AIC",
                                    ICCut = 2,
                                    data = speciesData,
                                    random = ~1|PlotID/uniTreeID, 
                                    control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  bestModel <- lme(fixed = as.formula(paste("logY~", paste(tempoutput$bestIDV, collapse = "+"))),
                   data = speciesData,
                   random = ~1|PlotID/uniTreeID, 
                   control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  
  cat("Species", indispecies, "is done. \n")
  if(indispecies == "All species"){
    modelSelectionSummaries <- tempoutput$modelSummary[, Species:=indispecies]
    bestFomula <- list(as.formula(paste("logY~", paste(tempoutput$bestIDV, collapse = "+"))))
    names(bestFomula) <- "All species"
    bestModels <- list(bestModel)
    names(bestModels) <- "All species"
  } else {
    bestFomula[[indispecies]] <- as.formula(paste("logY~", paste(tempoutput$bestIDV, collapse = "+")))
    modelSelectionSummaries <- rbind(modelSelectionSummaries, 
                                     tempoutput$modelSummary[, Species:=indispecies])
    bestModels[[indispecies]] <- bestModel
  }
}

  allFixedCoeff <- lapply(bestModels, function(x){
    data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                               conR2 = r.squaredGLMM(x)[2])]})
  




