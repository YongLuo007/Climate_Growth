rm(list = ls())
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
theplotidlebal <- unique(analysesData$PlotID)
thetreeidlebal <- unique(analysesData$uniTreeID)
analysesData[,':='(PlotID = factor(PlotID, levels = theplotidlebal, labels = theplotidlebal),
                   uniTreeID = factor(uniTreeID, levels = thetreeidlebal, 
                                      labels = thetreeidlebal))]

AllModels <- list()
AllResults <- data.table(Model = character(),
                         Formula = character(),
                         DIC = numeric(),
                         AIC = numeric(),
                         BIC = numeric(),
                         MarR2 = numeric(),
                         ConR2 = numeric(),
                         Species = character())

testSpecieses <- c("JP", "BS", "TA")
# testSpecieses <- c("PL", "AW", "SW", "SB", "PJ")
i <- 1
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
for(testspecies in testSpecieses){
  speciesData <- analysesData[Species == testspecies,]
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = Year-mean(Year),
                    logHctd = log(Hegyi)-mean(log(Hegyi)),
                    Dominancectd = Dominance_indiBiomass - mean(Dominance_indiBiomass))]
  
  modelselection <- mixedModelSelection(data = speciesData, DV = "logY", 
                                        maxInteraction = 3,
                                        IDV = c("logDBHctd", "Yearctd", "logHctd", "Dominancectd"),
                                        randomMethod = "both", slopeTerm = "Yearctd", 
                                        randomTerm = "PlotID/uniTreeID")
  
  AllResults <- rbind(AllResults, modelselection$modelSummary[, Species:=indispecies])
  names(modelselection$modelOutput) <- paste(indispecies, "_", 
                                             names(modelselection$modelOutput),
                                             sep = "")
  AllModels <- append(AllModels, modelselection$modelOutput)
  
}
rm(i, coeff, optim, testspecies, testSpecieses, theModel, analysesComponent,
   analysesData1, analysesDataOrg, DBHCut, minDBHdata, minFirstDBH, speciesData)

save.image(file.path(workPath, "Results","BiomassGR_mainResults.RData"))
