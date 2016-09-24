rm(list = ls())
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
theplotidlebal <- unique(analysesData$PlotID)
thetreeidlebal <- unique(analysesData$uniTreeID)
analysesData[,':='(PlotID = factor(PlotID, levels = theplotidlebal, labels = theplotidlebal),
                   uniTreeID = factor(uniTreeID, levels = thetreeidlebal, 
                                      labels = thetreeidlebal))]

OverallModels <- list()
OverallResults <- data.table(Species = character(),
                             MarR2 = numeric(),
                             ConR2 = numeric(),
                             Variable = character(),
                             Value = numeric(), Std.Error = numeric(),
                             DF = numeric(), t.value = numeric(),
                             p.value = numeric())
optim <- lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000)
testSpecieses <- c("JP", "BS", "TA")
# testSpecieses <- c("PL", "AW", "SW", "SB", "PJ")
i <- 1
for(testspecies in testSpecieses){
  speciesData <- analysesData[Species == testspecies,]
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = Year-mean(Year),
                    logHctd = log(Hegyi)-mean(log(Hegyi)),
                    Dominancectd = Dominance_indiBiomass - mean(Dominance_indiBiomass))]
  theModel <- lme(logY~logDBHctd+Yearctd+logHctd+Dominancectd+
                    logDBHctd:Yearctd+logDBHctd:logHctd+logDBHctd:Dominancectd+
                    Yearctd:logHctd+Yearctd:Dominancectd+logHctd:Dominancectd+
                    logDBHctd:Yearctd:logHctd+logDBHctd:Yearctd:Dominancectd+
                    Yearctd:logHctd:Dominancectd,
                  random = ~1+Yearctd|PlotID/uniTreeID,
                  data = speciesData,
                  control = optim)
  OverallModels[[i]] <- theModel
  names(OverallModels)[i] <- paste(testspecies, "_all", sep = "")
  i <- i+1
  coeff <- data.frame(summary(theModel)$tTable)
  coeff$Variable <- row.names(coeff)
  coeff <- data.table(coeff)[,':='(Species = testspecies, 
                                   MarR2 = as.numeric(r.squaredGLMM(theModel)[1]),
                                   ConR2 = as.numeric(r.squaredGLMM(theModel)[2]))]
  OverallResults <- rbind(OverallResults, coeff[,.(Species, MarR2, ConR2,
                                                   Variable,
                                                   Value, Std.Error, DF, t.value, p.value)])
}
rm(i, coeff, optim, testspecies, testSpecieses, theModel, analysesComponent,
   analysesData1, analysesDataOrg, DBHCut, minDBHdata, minFirstDBH, speciesData)

save.image(file.path(workPath, "Results","BiomassGR_mainResults.RData"))
