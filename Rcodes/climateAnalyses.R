rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table); library(parallel)
library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "newAllDataRescaledComp.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Other species")
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
bestClimateModels <- list()

climates <- c("ATA", "GSTA", "NONGSTA",
              "ACMIA", "GSCMIA", 
              "ACO2A")
for(indispecies in studySpecies){
  speciesData <- analysesData[DataType == indispecies,]
  for(indiclimate in climates){
    speciesData$climate <- c(speciesData[, indiclimate, with = FALSE])
    speciesData[,':='(logY = log(BiomassGR), 
                      logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                      Climatectd = climate-mean(climate),
                      logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                      logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                      logSActd = log(SA)-mean(log(SA)))]
    tempoutput <- mixedModelSelection(DV = "logY", 
                                      IDV = c("logDBHctd", "Climatectd", 
                                              "logIntraHctd", "logInterHctd",
                                              "logSActd"),
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
    
    if(indispecies == "All species" & indiclimate == "ATA"){
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
allFixedCoeff <- lapply(bestClimateModels, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})

save.image(file.path(workPath, "data", "ClimateModelSelection.RData"))
