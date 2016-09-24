rm(list = ls())
# this is the analyses for plot level biomass growth rate

workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
theplotidlebal <- unique(analysesData$PlotID)
thetreeidlebal <- unique(analysesData$uniTreeID)
analysesData[,':='(PlotID = factor(PlotID, levels = theplotidlebal, labels = theplotidlebal),
                   uniTreeID = factor(uniTreeID, levels = thetreeidlebal, 
                                      labels = thetreeidlebal))]

studyspecies <- c("JP", "TA", "BS")
FullModels <- list()
ReducedModels <- list()
AllResults <- data.table(Species = character(),
                         Model = character(),
                         DIC = numeric(),
                         MarR2 = numeric(),
                         ConR2 = numeric(),
                         Variable = character(),
                         Value = numeric(), Std.Error = numeric(),
                         DF = numeric(), t.value = numeric(),
                         p.value = numeric())
prevvari <- NULL
for(indispecies in studyspecies){
  
  speciesData <- analysesData[Species == indispecies,.(IniTotalBiomass = mean(IniBiomass),
                                                       FinTotalBiomass = mean(FinBiomass),
                                                       MeasureLength = mean(FinYear-IniYear),
                                                       meanH = mean(Hegyi),
                                                       meanDBH = mean(IniDBH),
                                                       meanDomi = mean(Dominance_indiBiomass)),
                              by = c("PlotID", "Year")]
  speciesData[, Plevel_BGR:=(FinTotalBiomass-IniTotalBiomass)/MeasureLength]
  # the full model
  speciesData[,':='(logY = log(Plevel_BGR), 
                    logDBHctd = log(meanDBH)-mean(log(meanDBH)), 
                    Yearctd = Year-mean(Year),
                    logHctd = log(meanH)-mean(log(meanH)),
                    Dominancectd = meanDomi - mean(meanDomi))]
  FullModel <- lme(logY~logDBHctd+Yearctd+logHctd+Dominancectd+
                     logDBHctd:Yearctd+logDBHctd:logHctd+logDBHctd:Dominancectd+
                     Yearctd:logHctd+Yearctd:Dominancectd+logHctd:Dominancectd+
                     logDBHctd:Yearctd:logHctd+logDBHctd:Yearctd:Dominancectd+
                     Yearctd:logHctd:Dominancectd,
                   random = ~1+Yearctd|PlotID,
                   data = speciesData)
  FullModels[[indispecies]] <- FullModel
  coeff <- data.frame(summary(FullModel)$tTable)
  coeff$Variable <- row.names(coeff)
  coeff <- data.table(coeff)[,':='(Species = indispecies, 
                                   Model = "Full",
                                   DIC = as.numeric(DIC(FullModel)),
                                   MarR2 = as.numeric(r.squaredGLMM(FullModel)[1]),
                                   ConR2 = as.numeric(r.squaredGLMM(FullModel)[2]))]
  AllResults <- rbind(AllResults, coeff[,.(Species, Model, DIC, MarR2, ConR2,
                                           Variable,
                                           Value, Std.Error, DF, t.value, p.value)])
  
  # reduced model
  for(i in 1:5){ # reduce model for 5 times
    variables <- coeff[p.value<0.05 & Variable != "(Intercept)", ]$Variable
    if(length(variables) != length(prevvari)){
      reducedFomu <- paste("logY~", variables[1], sep = "")
      if(length(variables)>=2){
        for(j in 2:length(variables)){
          reducedFomu <- paste(reducedFomu, "+", variables[j], sep = "")
        }
      }
      ReducedModel <- lme(fixed = as.formula(reducedFomu),
                          random = ~1+Yearctd|PlotID,
                          data = speciesData)
      ReducedModels[[indispecies]] <- ReducedModel
      coeff <- data.frame(summary(ReducedModel)$tTable)
      coeff$Variable <- row.names(coeff)
      coeff <- data.table(coeff)[,':='(Species = indispecies, 
                                       Model = paste("Reduced", i, sep = ""),
                                       DIC = as.numeric(DIC(ReducedModel)),
                                       MarR2 = as.numeric(r.squaredGLMM(ReducedModel)[1]),
                                       ConR2 = as.numeric(r.squaredGLMM(ReducedModel)[2]))]
      AllResults <- rbind(AllResults, coeff[,.(Species, Model, DIC, MarR2, ConR2,
                                               Variable,
                                               Value, Std.Error, DF, t.value, p.value)])
      prevvari <- variables
    }
  }
}

