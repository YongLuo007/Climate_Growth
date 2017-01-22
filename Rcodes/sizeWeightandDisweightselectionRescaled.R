rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
analysesData <- fread(file.path(workPath, "data", "newAllDataRescaledComp1.csv"))
analysesData[Species == "Other species", Species := "Minor species"]
analysesData[,TreeMT:=length(IniYear), by = c("Species", "uniTreeID")]
analysesData <- analysesData[TreeMT>=2, ]
analysesData <- analysesData[positiveGrowthTree == "yes",]

myFunction <- function(sizeWeight, disWeight, analysesData, workPath){
  output <- data.table(Species = character(), sizeWeight = character(),
                       disWeight = numeric(), HAIC = numeric(), IntraHAIC = numeric(),
                       InterHAIC = numeric())
  for(i in sizeWeight){
    for(j in disWeight){
      newCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
                                   paste("CompetitionData_DW",j, "_SW", i, ".csv", sep = "")))
      analysesDataAll <- setkey(analysesData, uniTreeID, IniYear)[setkey(newCIdata, uniTreeID, IniYear),
                                                                     nomatch = 0]
      
      for(indispecies in c("All species", "Jack pine", "Trembling aspen",
                           "Black spruce", "Minor species")){
        speciesData <- analysesDataAll[Species == indispecies,]
        speciesData[,':='(logY = log(BiomassGR),
                          logHctd = log(H+1)-mean(log(H+1)),
                          logIntraHctd = log(IntraH+1) - mean(log(IntraH+1)),
                          logInterHctd = log(InterH+1) - mean(log(InterH+1)))]
        HModel <- lme(logY~logHctd, random = ~1|uniTreeID,
                      data = speciesData,
                      control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
        IntraHModel <- lme(logY~logIntraHctd,
                       random = ~1|uniTreeID,
                       data = speciesData,
                       control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
        InterHModel <- lme(logY~logInterHctd,
                       random = ~1|uniTreeID,
                       data = speciesData,
                       control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
        aictable <- data.table(Species = indispecies, sizeWeight = i, 
                               disWeight = j, HAIC = AIC(HModel),
                               IntraHAIC = AIC(IntraHModel), 
                               InterHAIC = AIC(InterHModel))
        output <- rbind(output, aictable)
      }
    }
  }
  return(output)
}




inputWeights <- list()
m <- 1
for(i in seq(0, 10, by = 0.1)){
  for(j in seq(0, 2, by = 0.1)){
    inputWeights[[m]] <- c(i, j)
    m <- m+1
  }
}


analysesData[,':='(H = NULL, IntraH = NULL, InterH = NULL)]


cl <- parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl, c("lme", "AIC", "lmeControl", "data.table", "myFunction",
                              "setkey", "fread", "workPath", "analysesData"))

allresults <- parLapply(cl, inputWeights, 
          function(y) myFunction(sizeWeight = y[1], disWeight = y[2], 
                                 analysesData = analysesData, workPath = workPath))
stopCluster(cl)

for(i in 1:length(allresults)){
  if(i == 1){
    output <- allresults[[i]]
  } else {
    output <- rbind(output, allresults[[i]])
  }
}

write.csv(output, file.path(workPath, "Results", "bestAandBRandomUniTreeID_TwoCensusPositive.csv"), row.names = F)
# output <- fread(file.path(workPath, "Results", "newbestAandBRandomUniTreeID.csv"))
# output[,sizeWeight:=as.numeric(sizeWeight)]
# output <- output[sizeWeight<8.1,]
# output <- unique(output, by = c("Species", "sizeWeight", "disWeight"))
a <- melt(output, id.vars = c("Species", "sizeWeight", "disWeight"), 
          measure.vars = c("HAIC", "IntraHAIC", "InterHAIC"),
          value.name = "Value")

a[,minvalue:=min(Value), by = c("Species", "variable")]
bestWeightTable <- a[Value == minvalue,]

analysesData <- fread(file.path(workPath, "data", "newAllDataRescaledComp1.csv"))
studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce",
                  "Minor species")
analysesData[Species == "Other species", Species:="Minor species"]
analysesData[,':='(H = NULL, IntraH = NULL, InterH = NULL)]
for(indispecies in studySpecies){
  speciesData <- analysesData[Species == indispecies,]
  HsizeWeight <- bestWeightTable[Species == indispecies & variable == "HAIC",]$sizeWeight
  HdisWeight <- bestWeightTable[Species == indispecies & variable == "HAIC",]$disWeight
  HCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
                             paste("CompetitionData_DW",HdisWeight, "_SW", HsizeWeight, ".csv", sep = "")))
  
  IntraHsizeWeight <- bestWeightTable[Species == indispecies & variable == "IntraHAIC",]$sizeWeight
  IntraHdisWeight <- bestWeightTable[Species == indispecies & variable == "IntraHAIC",]$disWeight
  IntraHCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
                             paste("CompetitionData_DW",IntraHdisWeight, "_SW", IntraHsizeWeight, ".csv", sep = "")))
  
  InterHsizeWeight <- bestWeightTable[Species == indispecies & variable == "InterHAIC",]$sizeWeight
  InterHdisWeight <- bestWeightTable[Species == indispecies & variable == "InterHAIC",]$disWeight
  InterHCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
                                  paste("CompetitionData_DW",InterHdisWeight, "_SW", InterHsizeWeight, ".csv", sep = "")))
  allCIdata <- setkey(HCIdata[,.(uniTreeID, IniYear, H)], uniTreeID, IniYear)[setkey(IntraHCIdata[,.(uniTreeID, IniYear, IntraH)],
                                                                                     uniTreeID, IniYear), nomatch = 0]
  allCIdata <- setkey(allCIdata, uniTreeID, IniYear)[setkey(InterHCIdata[,.(uniTreeID, IniYear, InterH)],
                                                                                     uniTreeID, IniYear), nomatch = 0]
  
  speciesData <- setkey(speciesData, uniTreeID, IniYear)[setkey(allCIdata, uniTreeID, IniYear), 
                                                         nomatch = 0]
  if(indispecies == "All species"){
    newAllData <- speciesData
  } else {
    newAllData <- rbind(newAllData, speciesData)
  }
}

write.csv(newAllData, file.path(workPath, "data", "finalData.csv"), row.names = FALSE)





