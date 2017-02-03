rm(list=ls())
workPath <- "~/GitHub/Climate_Growth"

selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
load(file.path(workPath, "data", selectionMethod, "BestClimateModels.RData"))

allHFixedCoeff <- lapply(allHFixedCoeff, function(s) {
  s[, rn:=lapply(lapply(rn, function(x) sort(unlist(strsplit(x, ":", fixed = T)))),
                 function(z) paste(z, collapse = ":"))]})
a <- c("(Intercept)", "logDBHctd", "logSActd", "Climatectd", 
       "logHctd", "logDBHctd:logSActd",
       "Climatectd:logDBHctd", "logDBHctd:logHctd",
       "Climatectd:logSActd", 
       "logHctd:logSActd",
       "Climatectd:logHctd",
       "logInterHctd:logIntraHctd")
generralClimateCoeff <- data.table(rn = a[a %in% unique(unlist(lapply(allHFixedCoeff, function(s) s$rn)))])

ATA <- paste(studySpecies, "_ATA", sep = "")
GSTA <- paste(studySpecies, "_GSTA", sep = "")
NONGSTA <- paste(studySpecies, "_NONGSTA", sep = "")
ACMIA <- paste(studySpecies, "_ACMIA", sep = "")
GSCMIA <- paste(studySpecies, "_GSCMIA", sep = "")
ACO2A <- paste(studySpecies, "_ACO2A", sep = "")

for(i in c("ATA", "GSTA", "NONGSTA", "ACMIA",
           "GSCMIA", "ACO2A")){
  listnames <- get(i)
  k <- 1
  IndiclimateCoeff <- data.table::copy(generralClimateCoeff)
  for(j in listnames){
    indifixed <- allHFixedCoeff[j][[1]]
    indifixed <- indifixed[,.(rn = unlist(rn), Value = Value,
                              SE = `Std.Error`)]
    # if(i == "ACMIA" | i == "GSCMIA"){
    #   indifixed$Value[c(grep("Climatectd", indifixed$rn))] <- 10*indifixed$Value[c(grep("Climatectd", indifixed$rn))]
    #   indifixed$SE[c(grep("Climatectd", indifixed$rn))] <- 10*indifixed$SE[c(grep("Climatectd", indifixed$rn))]
    # }
    indifixed[,':='(Value = round(Value, 4), SE = round(SE, 4))]
    names(indifixed)[2:3] <- paste(studySpecies[k], names(indifixed)[2:3])
    k <- k+1
    IndiclimateCoeff <- left_join(IndiclimateCoeff, indifixed, by = "rn") %>% data.table
  }
  IndiclimateCoeff <- rbind(IndiclimateCoeff[1,][,rn:=i],
                            IndiclimateCoeff)
  
  IndiclimateCoeff[,rn:=gsub(pattern = "Climatectd", replacement = i, x = rn, fixed = TRUE)]
  if(i == "ATA"){
    alloutput <- IndiclimateCoeff
  } else {
    alloutput <- rbind(alloutput, IndiclimateCoeff)
  }
}
names(alloutput)[1] <- "Predictor"
alloutput[, Predictor:=gsub("ctd", "", Predictor)]
alloutput[Predictor == "(Intercept)", Predictor:="Intercept"]
alloutput[, Predictor:=gsub(":", " × ", Predictor, fixed = TRUE)]
workPath <- "~/GitHub/Climate_Growth"
write.csv(alloutput, file.path(workPath, "TablesFigures", 
                               "Table S3. fixedEffectForBestClimateModels_H.csv"),
          row.names = FALSE)

