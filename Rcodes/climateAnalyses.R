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
  for(indiclimate in climates){
    speciesData$climate <- c(speciesData[, indiclimate, with = FALSE])
    speciesData[,':='(logY = log(BiomassGR), 
                      logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                      Climatectd = climate-mean(climate),
                      logIntraHctd = log(IntraH1_3+1)-mean(log(IntraH1_3+1)),
                      logInterHctd = log(InterH0_4+1)-mean(log(InterH0_4+1)),
                      RBIctd = RBI - mean(RBI),
                      logSActd = log(SA)-mean(log(SA)),
                      logSBctd = log(PlotBiomass)-mean(log(PlotBiomass)))]
    tempoutput <- mixedModelSelection(DV = "logY", 
                                      IDV = c("logDBHctd", "Climatectd", "logIntraHctd", "logInterHctd",
                                              "RBIctd", "logSActd", "logSBctd"),
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
    bestClimateModels[[paste(indispecies, "_", indiclimate, sep = "")]] <- thebestmodel
    cat("Species", indispecies, "and climate", indiclimate, "is done. \n")
    
    if(indispecies == "All" & indiclimate == "ATA"){
      modelSelectionSummaries <- tempoutput$modelSummary[, ':='(Species = indispecies, Climate = indiclimate)]
      bestIDVs <- list(tempoutput$bestIDV)
      names(bestIDVs) <- paste(indispecies, "_", indiclimate, sep = "")
    } else {
      bestIDVs[[paste(indispecies, "_", indiclimate, sep = "")]] <- tempoutput$bestIDV
      modelSelectionSummaries <- rbind(modelSelectionSummaries, 
                                       tempoutput$modelSummary[, ':='(Species = indispecies, Climate = indiclimate)])
    }
  }
}

save.image(file.path(workPath, "data", "ClimateModelSelection.RData"))
