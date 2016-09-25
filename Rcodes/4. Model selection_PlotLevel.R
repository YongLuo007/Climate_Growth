rm(list = ls())
# this is the analyses for plot level biomass growth rate
library(dplyr); library(SpaDES); library(nlme); library(data.table)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
theplotidlebal <- unique(analysesData$PlotID)
thetreeidlebal <- unique(analysesData$uniTreeID)
analysesData[,':='(PlotID = factor(PlotID, levels = theplotidlebal, labels = theplotidlebal),
                   uniTreeID = factor(uniTreeID, levels = thetreeidlebal, 
                                      labels = thetreeidlebal))]

studyspecies <- c("JP", "TA", "BS")
AllModels <- list()
AllResults <- data.table(Model = character(),
                         IDV_Base = character(),
                         Direction = character(),
                         IDV_Processed = character(),
                         IDV_Length = numeric(),
                         All_Significant = logical(),
                         DIC = numeric(),
                         AIC = numeric(),
                         BIC = numeric(),
                         MarR2 = numeric(),
                         ConR2 = numeric(),
                         Species = character())

source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
for(indispecies in studyspecies){
  speciesData <- analysesData[Species == indispecies,.(IniTotalBiomass = mean(IniBiomass),
                                                       FinTotalBiomass = mean(FinBiomass),
                                                       MeasureLength = mean(FinYear-IniYear),
                                                       meanH = mean(Hegyi),
                                                       meanDBH = mean(IniDBH),
                                                       meanDomi = mean(Dominance_indiBiomass)),
                              by = c("PlotID", "Year")]
  speciesData[, Plevel_BGR:=(FinTotalBiomass-IniTotalBiomass)/MeasureLength]
  speciesData[,':='(logY = log(Plevel_BGR), 
                    logDBHctd = log(meanDBH)-mean(log(meanDBH)), 
                    Yearctd = Year-mean(Year),
                    logHctd = log(meanH)-mean(log(meanH)),
                    Dominancectd = log(meanDomi) - mean(log(meanDomi)))]
  
  
  modelselection <- mixedModelSelection(data = speciesData, 
                                        DV = "logY", 
                                        maxInteraction = 3,
                                        IDV = c("logDBHctd", "Yearctd", "logHctd", "Dominancectd"),
                                        random = ~1+Yearctd|PlotID, 
                                        control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000))
  
  AllResults <- rbind(AllResults, modelselection$modelSummary[, Species:=indispecies])
  names(modelselection$modelOutput) <- paste(indispecies, "_", 
                                             names(modelselection$modelOutput),
                                             sep = "")
  AllModels <- append(AllModels, modelselection$modelOutput)
}

for(indispecies in studyspecies){
  cat(indispecies, "\n")
  print(summary(AllModels[[paste(indispecies,"_ReducedModel1", sep = "")]])$tTable)
}


