rm(list = ls())
# summary results
library(dplyr); library(SpaDES); library(nlme); library(data.table); library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "theBestModels.RData"))
OverallResults <- data.table(Species = character(),
                             MarR2 = numeric(),
                             ConR2 = numeric(),
                             Variable = character(),
                             Value = numeric(), Std.Error = numeric(),
                             DF = numeric(), t.value = numeric(),
                             p.value = numeric())
for(indispecies in c("JP", "TA", "BS")){
  themodel <- thebestmodel[[paste(indispecies, "_bestModel", sep = "")]]
  
  coeff <- data.frame(summary(themodel)$tTable)
  coeff$Variable <- row.names(coeff)
  coeff <- data.table(coeff)[,':='(Species = indispecies,
                                   MarR2 = as.numeric(r.squaredGLMM(themodel)[1]),
                                   ConR2 = as.numeric(r.squaredGLMM(themodel)[2]))]
  OverallResults <- rbind(OverallResults, coeff[,.(Species, MarR2, ConR2, 
                                                   Variable,
                                                   Value, Std.Error, DF, t.value, p.value)])
}

write.csv(OverallResults, file.path(workPath, "Results", "Year_model_results.csv"),
          row.names = FALSE)

