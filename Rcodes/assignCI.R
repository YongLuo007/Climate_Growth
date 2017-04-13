rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}

firstRun <- TRUE
if(firstRun){
  analysesDataOrg <- fread(file.path(workPath, "data", "Year10Analyses", 
                                  "MBdatafinal.csv"))
  analysesData <- rbind(data.table::copy(analysesDataOrg)[,Species_Group := "All species"],
                        data.table::copy(analysesDataOrg))
  analysesData <- analysesData[allCensusLiveTree == "yes",]
  myFunction <- function(sizeWeight, disWeight, analysesData, workPath){
    output <- data.table(Species = character(), sizeWeight = character(),
                         disWeight = numeric(), HAIC = numeric())
    for(i in sizeWeight){
      for(j in disWeight){
        newCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
                                     paste("CompetitionData_DW",j, "_SW", i, ".csv", sep = "")))
        setnames(newCIdata, "IniYear", "MidYear")
        analysesDataAll <- setkey(analysesData, uniTreeID, MidYear)[setkey(newCIdata, uniTreeID, MidYear),
                                                                    nomatch = 0]
        
        for(indispecies in c("All species", "Jack pine", "Trembling aspen",
                             "Black spruce", "Minor species")){
          speciesData <- analysesDataAll[Species_Group == indispecies,]
          minABGR <- abs(round(min(speciesData$BiomassGR), 3))+0.01
          speciesData[,':='(logY = log(BiomassGR+minABGR),
                            logHctd = log(H)-mean(log(H)))]
          HModel <- lme(logY~logHctd,
                        random = ~1|uniTreeID,
                        data = speciesData,
                        control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
          aictable <- data.table(Species = indispecies, sizeWeight = i, 
                                 disWeight = j, HAIC = AIC(HModel))
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
  
  
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c("lme", "AIC", "lmeControl", "data.table", "myFunction",
                                "setkey", "fread", "workPath", "analysesData", "setnames"))
  
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
  
  write.csv(output, file.path(workPath, "data", "Year10Analyses",
                              paste("bestAandB.csv", sep = "")),
            row.names = F)
} else {
  
  output <- fread(file.path(file.path(workPath, "data", "Year10Analyses",
                                      paste("bestAandB.csv", sep = ""))))
}


a <- melt(output, id.vars = c("Species", "sizeWeight", "disWeight"), 
          measure.vars = c("HAIC"),
          value.name = "Value")

a[,minvalue:=min(Value), by = c("Species", "variable")]
bestWeightTable <- a[Value == minvalue,]

analysesDataOrg <- fread(file.path(workPath, "data", "Year10Analyses", 
                                   "MBdatafinal.csv"))
analysesData <- rbind(data.table::copy(analysesDataOrg)[,Species_Group := "All species"],
                      data.table::copy(analysesDataOrg))

studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce",
                  "Minor species")


for(indispecies in studySpecies){
  speciesData <- analysesData[Species_Group == indispecies,]
  HsizeWeight <- bestWeightTable[Species == indispecies & variable == "HAIC",]$sizeWeight
  HdisWeight <- bestWeightTable[Species == indispecies & variable == "HAIC",]$disWeight
  HCIdata <- fread(file.path(workPath, "data", "AllCompetitionData",
                             paste("CompetitionData_DW",HdisWeight, "_SW", HsizeWeight, ".csv", sep = "")))
  setnames(HCIdata, "IniYear", "MidYear")
  speciesData <- setkey(speciesData, uniTreeID, MidYear)[setkey(HCIdata, uniTreeID, MidYear), 
                                                         nomatch = 0]
  speciesData[, c("IntraH", "InterH"):=NULL]
  
  if(indispecies == "All species"){
    newAllData <- speciesData
  } else {
    newAllData <- rbind(newAllData, speciesData)
  }
}

write.csv(newAllData, file.path(workPath, "data", "Year10Analyses", "finalData10.csv"),
          row.names = FALSE)
