rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
studySpecies <- c("All", "JP", "BS", "TA", "Other")
majorSpecies <- c("JP", "TA", "BS")
analysesData[,':='(IntraH = Hegyi*IntraHegyiRatio, InterH = Hegyi*(1-IntraHegyiRatio))]


climates <- c("ATA", "GSTA", "NONGSTA",
              "ACMIA", "GSCMIA", "NONGSCMIA",
              "ACO2A", "GSCO2A", "NONGSCO2A")
allClimateModels <- list()
for(indispecies in studySpecies){
  if(indispecies == "All"){
    speciesData <- data.table::copy(analysesData)
  } else if(indispecies == "Other"){
    speciesData <- analysesData[!(Species %in% majorSpecies),]
  } else {
    speciesData <- analysesData[Species == indispecies,]
  }
  for(indiclimate in climates){
    if(indiclimate != "ATA"){ set(speciesData, ,"climate", NULL)}
    setnames(speciesData, indiclimate, "climate")
    speciesData[,':='(logY = log(BiomassGR), 
                      logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                      climatectd = climate-mean(climate),
                      logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                      logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                      RBIctd = RBI - mean(RBI))]
    fullthemodel <- lme(logY~logDBHctd+climatectd+logIntraHctd+logInterHctd+RBIctd+
                      climatectd:logDBHctd+logDBHctd:logIntraHctd+
                      logDBHctd:logInterHctd+logDBHctd:RBIctd+
                      climatectd:logIntraHctd+climatectd:logInterHctd+
                      climatectd:RBIctd+logInterHctd:logIntraHctd+
                      logIntraHctd:RBIctd+logInterHctd:RBIctd+
                      climatectd:logDBHctd:logIntraHctd+
                      climatectd:logDBHctd:logInterHctd+
                      climatectd:logDBHctd:RBIctd+
                      logDBHctd:logInterHctd:logIntraHctd+
                      logDBHctd:logIntraHctd:RBIctd+
                      logDBHctd:logInterHctd:RBIctd+
                      climatectd:logInterHctd:logIntraHctd+
                      climatectd:logIntraHctd:RBIctd+
                      climatectd:logInterHctd:RBIctd+
                      logInterHctd:logIntraHctd:RBIctd,
                    random = ~1|PlotID/uniTreeID, 
                    control = lmeControl(opt="optim", maxIter=20000, msMaxIter = 20000),
                    data = speciesData)
    signIDV <- data.table(summary(fullthemodel)$tTable, keep.rownames = T)[`p-value` < 0.05 & rn != "(Intercept)",]$rn
    reducedFormu <- as.formula(paste("logY~", paste(signIDV, collapse = "+")))
    reducedModel <- lme(fixed = reducedFormu,
                        random = ~1|PlotID/uniTreeID, 
                        control = lmeControl(opt="optim", maxIter=20000, msMaxIter = 20000),
                        data = speciesData)
    newSigNIDV <- data.table(summary(reducedModel)$tTable, keep.rownames = T)[`p-value` < 0.05 & rn != "(Intercept)",]$rn
    for(i in 1:5){
      if(length(signIDV) != length(newSigNIDV)){
        signIDV <- newSigNIDV
        reducedFormu <- as.formula(paste("logY~", paste(newSigNIDV, collapse = "+")))
        reducedModel <- lme(fixed = reducedFormu,
                            random = ~1|PlotID/uniTreeID, 
                            control = lmeControl(opt="optim", maxIter=20000, msMaxIter = 20000),
                            data = speciesData)
        newSigNIDV <- data.table(summary(reducedModel)$tTable, keep.rownames = T)[`p-value` < 0.05 & rn != "(Intercept)",]$rn
      }
    }
    allClimateModels[[paste(indispecies, "_", indiclimate, sep = "")]] <- reducedModel
  }
}

save.image(file.path(workPath, "data", "ClimateResults.RData"))
