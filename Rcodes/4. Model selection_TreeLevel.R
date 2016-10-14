rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
AllResults <- data.table(Model = character(),
                         Formula = character(),
                         AIC = numeric(),
                         deltaAIC = numeric(),
                         MarR2 = numeric(),
                         ConR2 = numeric(),
                         Species = character())
analysesData[,':='(IntraH = Hegyi*IntraHegyiRatio, InterH = Hegyi*(1-IntraHegyiRatio))]

studySpecies <- c("All", "JP", "TA", "BS", "Other")[1]
majorSpecies <- c("JP", "TA", "BS")
i <- 1
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
                    Yearctd = Year-mean(Year),
                    logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                    logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                    RBIctd = RBI - mean(RBI))]
  
  modelselection <- mixedModelSelection(DV = "logY", 
                                        IDV = c("logDBHctd", "Yearctd", "logIntraHctd",
                                                "logInterHctd", "RBIctd"),
                                        maxInteraction = 3,
                                        ICTerm = "AIC",
                                        data = speciesData, 
                                        random = ~1+Yearctd|PlotID/uniTreeID, 
                                        control = lmeControl(opt="optim", maxIter=20000,
                                                             msMaxIter = 20000))
  
  AllResults <- rbind(AllResults, modelselection$modelSummary[, Species:=indispecies])
}

save.image(file.path(workPath, "data",
                     paste(indispecies, "_BiomassGR_Year_modelselection.RData", sep = "")))
