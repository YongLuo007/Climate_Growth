rm(list = ls())
# summary results
library(dplyr); library(SpaDES); library(nlme); library(data.table); library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"

# load(file.path(workPath, "data", "BiomassGR_Year_Model_Final.RData"))
load(file.path(workPath, "data", "theBestModels.RData"))
allmodels <- allbestmodels
OverallResults <- data.table(Species = character(),
                             MarR2 = numeric(),
                             ConR2 = numeric(),
                             Variable = character(),
                             Value = numeric(), Std.Error = numeric(),
                             DF = numeric(), t_value = numeric(),
                             p_value = numeric())
for(indispecies in c("All", "JP", "TA", "BS", "Other")){
  themodel <- allmodels[[indispecies]]
  
  coeff <- data.table(summary(themodel)$tTable, keep.rownames = TRUE)
  names(coeff)[5:6] <- c("t_value", "p_value")
  coeff <- data.table(coeff)[,':='(Species = indispecies,
                                   MarR2 = as.numeric(r.squaredGLMM(themodel)[1]),
                                   ConR2 = as.numeric(r.squaredGLMM(themodel)[2]))]
  OverallResults <- rbind(OverallResults, coeff[,.(Species, MarR2, ConR2, 
                                                   Variable = rn,
                                                   Value, Std.Error, DF, t_value, p_value)])
}

write.csv(OverallResults, file.path(workPath, "Results", "Year_model_results.csv"),
          row.names = FALSE)

