rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(lme4)
if(as.character(Sys.info()[1]) == "Windows"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}

selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
load(file.path(workPath, "data", selectionMethod, "BestYearModels.RData"))

selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
newFormula <- lapply(indiHbestFormulas, function(s) paste(as.character(s)[2], "~",
                                                         as.character(s)[3], "+(1|PlotID/uniTreeID)",
                                                         sep = ""))
allCoeffs <- data.frame(Variable = c("Intercept", "logDBH", "logSA", "Year", 
                                     "logIntraH", "logInterH", "logDBH:logSA",
                                     "logDBH:Year", "logDBH:logIntraH",
                                     "logDBH:logInterH", "logSA:Year", 
                                     "logIntraH:logSA", "logInterH:logSA",
                                     "logIntraH:Year", "logInterH:Year",
                                     "logInterH:logIntraH"),
                        stringsAsFactors = FALSE)

for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies,]
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = Year-mean(Year),
                    logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                    logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                    logSActd = log(IniFA+2.5)-mean(log(IniFA+2.5)))]
  allHbestModel <- lmer(formula = as.formula(newFormula[[indispecies]]),
                        data = speciesData)
  indiANOVA <- as.data.table(anova(allHbestModel),
                             keep.rownames = TRUE)[,rn:=gsub("ctd", "", rn)]
  indiANOVA <- indiANOVA[,.(Variable = rn, SS = round(`Sum Sq`, 2),
                            F = round(`F value`, 2))]
  indiANOVA[, Variable:=unlist(lapply(Variable,
                                      function(s) 
                                        paste(sort(unlist(strsplit(s, ":", fixed = TRUE))), collapse = ":")))]
  names(indiANOVA)[2:3] <- paste(indispecies, "_", names(indiANOVA)[2:3], sep = "")
  indiCoeff <- as.data.table(summary(allHbestModel)$coefficients,
                             keep.rownames = TRUE,
                             stringsAsFactor = FALSE)[,rn:=gsub("ctd", "", rn)]
  indiCoeff[rn=="(Intercept)", rn:="Intercept"]
  indiCoeff <- indiCoeff[,.(Variable = rn, 
                            Value = paste(round(Estimate, 4),
                                          "(", round(`Std. Error`, 4),
                                          ")", sep = ""))]
  indiCoeff[, Variable:=unlist(lapply(Variable,
                                      function(s) 
                                        paste(sort(unlist(strsplit(s, ":", fixed = TRUE))), collapse = ":")))]
  names(indiCoeff)[2] <- paste(indispecies, "_", names(indiCoeff)[2], sep = "")
  allCoeffs <- dplyr::left_join(allCoeffs, indiCoeff, by = "Variable")
  allCoeffs <- dplyr::left_join(allCoeffs, indiANOVA, by = "Variable")
}

workPath <- "~/Github/Climate_Growth"
write.csv(allCoeffs, file.path(workPath, "TablesFigures", "Table S2. fixedEffect_BestYearModels_IntraInterH.csv"),
          row.names = FALSE)
