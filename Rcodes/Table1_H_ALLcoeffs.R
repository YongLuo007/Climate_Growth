rm(list = ls())
library(dplyr); library(SpaDES); library(lme4); library(data.table);library(MuMIn)
library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
selectionMethod <- "Year10Analyses"
analysesDataOrg <- fread(file.path(workPath, "data", selectionMethod, "finalData10.csv"))
analysesDataAll <- analysesDataOrg[allCensusLiveTree == "yes",]
allCoeffs <- data.frame(Variable = c("Intercept", "logDBH", "logSA", "logH", "Year",
                                     "logDBH:logSA", "logDBH:logH", "logDBH:Year",
                                     "logH:logSA", "logSA:Year",
                                     "logH:Year"),
                        stringsAsFactors = FALSE)
studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Minor species")

for(indispecies in studySpecies){
  speciesDataAll <- analysesDataAll[Species_Group == indispecies,]
  minABGR <- round(abs(min(speciesDataAll$BiomassGR)), 3)+0.01
  speciesDataAll[,':='(logY = log(BiomassGR+minABGR), 
                       logDBHctd = log(MidDBH)-mean(log(MidDBH)), 
                       Yearctd = MidYear-mean(MidYear),
                       logHctd = log(H)-mean(log(H)),
                       logSActd = log(MidFA)-mean(log(MidFA)))]
  
  fullModelAll <- lmer(logY~logDBHctd+Yearctd+logHctd+logSActd+
                         logDBHctd:Yearctd+logDBHctd:logHctd+logDBHctd:logSActd+
                         Yearctd:logHctd+Yearctd:logSActd+
                         logHctd:logSActd+(Yearctd+1|PlotID/uniTreeID),
                       data = speciesDataAll)
  
  
  indiANOVA <- as.data.table(anova(fullModelAll),
                             keep.rownames = TRUE)[,rn:=gsub("ctd", "", rn)]
  indiANOVA <- indiANOVA[,.(Variable = rn, SS = round(`Sum Sq`, 2))]
  indiANOVA[, Variable:=unlist(lapply(Variable,
                                      function(s) 
                                        paste(sort(unlist(strsplit(s, ":", fixed = TRUE))),
                                              collapse = ":")))]
  names(indiANOVA)[2] <- paste(indispecies, "_", names(indiANOVA)[2], sep = "")
  indiCoeff <- as.data.table(summary(fullModelAll)$coefficients,
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
  cat("Species", indispecies, "is done. \n")
  rm(fullModelAll)
}



workPath <- "~/Github/Climate_Growth"
write.csv(allCoeffs, file.path(workPath, "TablesFigures", "Table1.csv"),
          row.names = FALSE)
