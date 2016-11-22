rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table); library(parallel)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
studySpecies <- c("All", "JP", "BS", "TA", "Other")
majorSpecies <- c("JP", "TA", "BS")
bestClimateModels <- list()

climates <- c("ATA", "GSTA", "NONGSTA",
              "ACMIA", "GSCMIA", 
              "ACO2A")
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))

for(indispecies in studySpecies){
  if(indispecies == "All"){
    speciesData <- data.table::copy(analysesData)
  } else if(indispecies == "Other"){
    speciesData <- analysesData[!(Species %in% majorSpecies),]
  } else {
    speciesData <- analysesData[Species == indispecies,]
  }
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    ATActd = ATA-mean(ATA),
                    ACMIActd = ACMIA-mean(ACMIA),
                    ACO2Actd = ACO2A-mean(ACO2A),
                    logIntraHctd = log(IntraH1_3+1)-mean(log(IntraH1_3+1)),
                    logInterHctd = log(InterH0_4+1)-mean(log(InterH0_4+1)),
                    RBIctd = RBI - mean(RBI),
                    logSActd = log(SA)-mean(log(SA)),
                    logSBctd = log(PlotBiomass)-mean(log(PlotBiomass)))]
  tempoutput <- mixedModelSelection(DV = "logY", 
                                    IDV = c("logDBHctd", "ATActd", "ACMIActd", "ACO2Actd", "logIntraHctd",
                                            "logInterHctd", "RBIctd", "logSActd", "logSBctd"),
                                    maxInteraction = 2,
                                    ICTerm = "AIC",
                                    ICCut = 2,
                                    data = speciesData,
                                    random = ~1|PlotID/uniTreeID, 
                                    control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  thebestmodel <- lme(as.formula(paste("logY~", paste(tempoutput$bestIDV, collapse = "+"))),
                      data = speciesData,
                      random = ~1|PlotID/uniTreeID, 
                      control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  bestClimateModels[[indispecies]] <- thebestmodel
  cat("Species", indispecies, "is done. \n")
  
  if(indispecies == "All"){
    modelSelectionSummaries <- tempoutput$modelSummary[, ':='(Species = indispecies)]
    bestIDVs <- list(tempoutput$bestIDV)
    names(bestIDVs) <- indispecies
  } else {
    bestIDVs[[indispecies]] <- tempoutput$bestIDV
    modelSelectionSummaries <- rbind(modelSelectionSummaries, 
                                     tempoutput$modelSummary[, ':='(Species = indispecies)])
  }
}

save.image(file.path(workPath, "data", "ClimateModelTogetherSelection.RData"))
