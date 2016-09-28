rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
theplotidlebal <- unique(analysesData$PlotID)
thetreeidlebal <- unique(analysesData$uniTreeID)
analysesData[,':='(PlotID = factor(PlotID, levels = theplotidlebal, labels = theplotidlebal),
                   uniTreeID = factor(uniTreeID, levels = thetreeidlebal, 
                                      labels = thetreeidlebal))]


AllResults <- data.table(Model = character(),
                         Formula = character(),
                         AIC = numeric(),
                         deltaAIC = numeric(),
                         MarR2 = numeric(),
                         ConR2 = numeric(),
                         Species = character())


testSpecieses <- c("JP", "BS", "TA")
# testSpecieses <- c("PL", "AW", "SW", "SB", "PJ")
i <- 1
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
for(indispecies in testSpecieses[1]){
  speciesData <- analysesData[Species == indispecies,]
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = Year-mean(Year),
                    logHctd = log(Hegyi)-mean(log(Hegyi)),
                    RBIctd = RBI - mean(RBI))]
  
  modelselection <- mixedModelSelection(DV = "logY", 
                                        IDV = c("logDBHctd", "Yearctd", "logHctd", "RBIctd"),
                                        maxInteraction = 3,
                                        ICTerm = "AIC",
                                        data = speciesData, 
                                        random = ~1+Yearctd|PlotID/uniTreeID, 
                                        control = lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000))
  
  AllResults <- rbind(AllResults, modelselection$modelSummary[, Species:=indispecies])
}

save.image(file.path(workPath, "Results",
                     paste(indispecies, "_BiomassGR_Year_modelselection.RData", sep = "")))
