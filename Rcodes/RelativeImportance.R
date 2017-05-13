rm(list = ls())
library(relaimpo);library(data.table);library(ggplot2)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
selectionMethod <- "Year10Analyses"
analysesDataOrg <- fread(file.path(workPath, "data", selectionMethod, "finalData10.csv"))
analysesData <- analysesDataOrg[allCensusLiveTree == "yes",]

studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Minor species")

allIDVSallH <- rbind(data.table(Variable = c("logDBH", "logSA", "logH", "logDBHlogSA",
                                             "logHlogSA", "logDBHlogH"),
                                Group = "Ontogeny+Competition"),
                     data.table(Variable = c("Year"),
                                Group = "Climate"),
                     data.table(Variable = c("logDBHYear","logSAYear"),
                                Group = "Ontogeny+Climate"),
                     data.table(Variable = c("logHYear"),
                                Group = "Competition+Climate"))
groupnames <- unique(allIDVSallH$Group)
grouplistallH <- list()
for(i in groupnames){
  grouplistallH[[i]] <- allIDVSallH[Group == i,]$Variable
}



for(indispecies in studySpecies){
  speciesData <- analysesData[Species_Group == indispecies, ]
  minABGR <- round(abs(min(speciesData$BiomassGR)), 3)+0.01
  speciesData[,':='(logY = log(BiomassGR+minABGR), 
                    logDBH = log(MidDBH)-mean(log(MidDBH)), 
                    Year = MidYear-mean(MidYear),
                    logH = log(H)-mean(log(H)),
                    logSA = log(MidFA)-mean(log(MidFA)))]
  
  speciesData[,':='(logDBHlogH = logDBH*logH,
                    logDBHlogSA = logDBH*logSA, 
                    logDBHYear = logDBH*Year, 
                    logHlogSA = logH*logSA,
                    logHYear = logH*Year,
                    logSAYear = logSA*Year)]
  suppressWarnings(des <- svydesign(id=~PlotID/uniTreeID,data = speciesData))
  bt1 <- boot.relimp(logY~logDBH+Year+logH+logSA+
                       logDBHYear+logDBHlogH+logDBHlogSA+
                       logHYear+logSAYear+
                       logHlogSA,
                     data = speciesData,
                     type = c("lmg"),
                     b = 1000,
                     design = des,
                     # groups = grouplistallH,
                     # groupnames = groupnames,
                     rela = TRUE)
  bt <- booteval.relimp(bt1, level = c(0.95))
  tempfixedLmg <- attributes(bt)
  temptable <- data.table(Species = indispecies, Group = tempfixedLmg$namen[-1],
                          importance = tempfixedLmg$lmg, 
                          importanceLower = as.numeric(tempfixedLmg$lmg.lower),
                          importanceUpper = as.numeric(tempfixedLmg$lmg.upper))
  cat("Species: ", indispecies, "is done. \n")
  if(indispecies == "All species"){
    importanceTable <- temptable
  } else {
    importanceTable <- rbind(importanceTable, temptable)
  }
  cat("Species", indispecies, "is done. \n")
}
importanceTable[,':='(importance=importance*100, 
                      importanceLower=importanceLower*100,
                      importanceUpper=importanceUpper*100)]
write.csv(importanceTable, file.path(workPath, "TablesFigures", "importanceTable.csv"), row.names = F)


