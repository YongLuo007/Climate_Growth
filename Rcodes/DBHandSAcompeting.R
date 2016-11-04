rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table); library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
studySpecies <- c("All", "JP", "BS", "TA", "Other")
AICtable <- data.table(Species = character(), AICDBH = numeric(),
                       AICSA = numeric(), AICDBHwithSA = numeric(),
                       AICH = numeric(),
                       AICSB = numeric(), AICHwithSB = numeric(),
                       AICHwithSA = numeric(), AICHwithSASB = numeric())
for(indispecies in studySpecies){
  if(indispecies == "All"){
    speciesData <- analysesData
  } else if(indispecies == "Other"){
    speciesData <- analysesData[!(Species %in% c("JP", "BS", "TA")),]
  } else {
    speciesData <- analysesData[Species == indispecies,]
  }
  speciesData[,':='(logY = log(BiomassGR), 
                    logSActd = log(SA)-mean(log(SA)),
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)),
                    logIntraHctd = log(IntraH1_3+1)-mean(log(IntraH1_3+1)),
                    logInterHctd = log(InterH0_4+1)-mean(log(InterH0_4+1)),
                    logSBctd = log(PlotBiomass)-mean(log(PlotBiomass)))]
  DBHmodel <- lme(logY~logDBHctd,
                  random = ~1|PlotID,
                  data = speciesData,
                  control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  Hmodel <- lme(logY~logIntraHctd+logInterHctd+logIntraHctd:logInterHctd,
                random = ~1|PlotID,
                data = speciesData,
                control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  SAmodel <- lme(logY~logSActd,
                 random = ~1|PlotID,
                 data = speciesData,
                 control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  SBmodel <- lme(logY~logSBctd,
                 random = ~1|PlotID,
                 data = speciesData,
                 control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  DBHandSAmodel <- lme(logY~logDBHctd+logSActd,
                       random = ~1|PlotID,
                       data = speciesData,
                       control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  HandSBmodel <- lme(logY~logIntraHctd+logInterHctd+logIntraHctd:logInterHctd + logSBctd,
                     random = ~1|PlotID,
                     data = speciesData,
                     control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  HandSAmodel <- lme(logY~logIntraHctd+logInterHctd+logIntraHctd:logInterHctd + logSActd,
                     random = ~1|PlotID,
                     data = speciesData,
                     control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  HandSASBmodel <- lme(logY~logIntraHctd+logInterHctd+logIntraHctd:logInterHctd + logSBctd+
                       logSActd,
                     random = ~1|PlotID,
                     data = speciesData,
                     control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  AICtable <- rbind(AICtable,
                    data.table(Species = indispecies, 
                               AICDBH = AIC(DBHmodel), 
                               AICSA = AIC(SAmodel), 
                               AICDBHwithSA = AIC(DBHandSAmodel), 
                               AICH = AIC(Hmodel), 
                               AICSB = AIC(SBmodel), 
                               AICHwithSB = AIC(HandSBmodel),
                               AICHwithSA = AIC(HandSAmodel),
                               AICHwithSASB = AIC(HandSASBmodel)))
}



