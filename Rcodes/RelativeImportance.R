rm(list = ls())
library(relaimpo);library(data.table);library(ggplot2)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
load(file.path(workPath, "data", "finalYearModels168plotsAllCensusPositiveTrees.RData"))
bestallHIDVs <- lapply(allHbestFormulas, function(s) as.character(s)[[3]])
bestallHIDVs <- lapply(bestallHIDVs, function(x) {unlist(strsplit(x, " + ", fixed = TRUE))})
bestallHIDVsnew <- list()
for(i in 1:length(bestallHIDVs)){
  indibestIDV <- bestallHIDVs[[i]]
  indibestIDV <- unlist(lapply(lapply(lapply(lapply(indibestIDV, function(x) strsplit(x, ":")), function(s) unlist(s)),
                                      function(f) sort(f)), function(d) paste(d, collapse = "")))
  bestallHIDVsnew[[i]] <- indibestIDV
}
bestallHIDVsnew <- lapply(bestallHIDVsnew, function(x) gsub("ctd", "", x))
names(bestallHIDVsnew) <- names(allHbestFormulas)
allHformulas <- lapply(bestallHIDVsnew, function(x) paste("logY~", paste(x, collapse = "+")))
names(allHformulas) <- names(allHbestFormulas)

bestindiHIDVs <- lapply(indiHbestFormulas, function(s) as.character(s)[[3]])
bestindiHIDVs <- lapply(bestindiHIDVs, function(x) {unlist(strsplit(x, " + ", fixed = TRUE))})
bestindiHIDVsnew <- list()
for(i in 1:length(bestindiHIDVs)){
  indibestIDV <- bestindiHIDVs[[i]]
  indibestIDV <- unlist(lapply(lapply(lapply(lapply(indibestIDV, function(x) strsplit(x, ":")), function(s) unlist(s)),
                                      function(f) sort(f)), function(d) paste(d, collapse = "")))
  bestindiHIDVsnew[[i]] <- indibestIDV
}
bestindiHIDVsnew <- lapply(bestindiHIDVsnew, function(x) gsub("ctd", "", x))
names(bestindiHIDVsnew) <- names(indiHbestFormulas)
indiHformulas <- lapply(bestindiHIDVsnew, function(x) paste("logY~", paste(x, collapse = "+")))
names(indiHformulas) <- names(indiHbestFormulas)


ontogeny <- c("logDBH", "logSA", "logDBHlogSA")
competition <- c("logIntraH", "logInterH", "logInterHlogIntraH", "logH")
climate <- "Year"
ontogenyCompetition <- c("logDBHlogInterH", "logIntraHlogSA", 
                         "logInterHlogSA", "logDBHlogIntraH",
                         "logHlogSA", "logDBHlogH")
ontogenyClimate <- c("logDBHYear","logSAYear")
competitionClimate <- c("logIntraHYear", 
                        "logInterHYear",
                        "logHYear")

for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies, ][
    ,':='(logY = log(BiomassGR),
          logDBH = log(IniDBH)-mean(log(IniDBH)),
          logSA = log(IniFA+2.5) - mean(log(IniFA+2.5)),
          logH = log(H+1)-mean(log(H+1)),
          logIntraH = log(IntraH+1) - mean(log(IntraH+1)),
          logInterH = log(InterH+1)- mean(log(InterH+1)),
          Year = Year-mean(Year))]
  speciesData[,':='(logDBHlogInterH = logDBH*logInterH, 
                    logDBHlogIntraH = logDBH*logIntraH,
                    logDBHlogH = logDBH*logH,
                    logDBHlogSA = logDBH*logSA, 
                    logDBHYear = logDBH*Year, 
                    logInterHlogIntraH = logInterH*logIntraH, 
                    logInterHlogSA = logInterH*logSA,
                    logInterHYear = logInterH*Year,
                    logHlogSA = logH*logSA,
                    logIntraHlogSA = logIntraH*logSA, 
                    logIntraHYear = logIntraH*Year,
                    logHYear = logH*Year,
                    logSAYear = logSA*Year)]
  allIDVSallH <- data.table(Variable = bestallHIDVsnew[[indispecies]])
  allIDVSallH[Variable %in% ontogeny, Group:="Ontogeny"]
  allIDVSallH[Variable %in% competition, Group:="Competition"]
  allIDVSallH[Variable %in% climate, Group:="Climate"]
  allIDVSallH[Variable %in% ontogenyCompetition, Group:="Ontogeny+Competition"]
  allIDVSallH[Variable %in% ontogenyClimate, Group:="Ontogeny+Climate"]
  allIDVSallH[Variable %in% competitionClimate, Group:="Competition+Climate"]
  groupnamesallH <- unique(allIDVSallH$Group)
  grouplistallH <- list()
  for(i in groupnames){
    grouplistallH[[i]] <- allIDVSallH[Group == i,]$Variable
  }
  suppressWarnings(des <- svydesign(id=~uniTreeID,data = speciesData))
  bt1 <- boot.relimp(as.formula(allHformulas[[indispecies]]),
                     data = speciesData,
                     type = c("lmg"),
                     b = 2000,
                     design = des,
                     groups = grouplistallH,
                     groupnames = groupnames,
                     rela = TRUE)
  bt <- booteval.relimp(bt1, level = c(0.95))
  tempfixedLmg <- attributes(bt)
  temptable <- data.table(Species = indispecies, Group = tempfixedLmg$namen[-1],
                          importance = tempfixedLmg$lmg, 
                          importanceLower = as.numeric(tempfixedLmg$lmg.lower),
                          importanceUpper = as.numeric(tempfixedLmg$lmg.upper))
  temptable[Group %in% ontogeny, Group:="Ontogeny"]
  temptable[Group %in% competition, Group:="Competition"]
  temptable[Group %in% climate, Group:="Climate"]
  temptable[Group %in% ontogenyCompetition, Group:="Ontogeny+Competition"]
  temptable[Group %in% ontogenyClimate, Group:="Ontogeny+Climate"]
  temptable[Group %in% competitionClimate, Group:="Competition+Climate"]
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
write.csv(importanceTable, file.path(workPath, "data", "importanceTable.csv"), row.names = F)


