rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra);library(grid)
library(raster);library(maptools);library(rgeos)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
load(file.path(workPath, "data",selectionMethod,
               "BestYearModels.RData"))
workPath <- "~/GitHub/Climate_Growth"

analysesData <- fread(file.path(workPath, "data", selectionMethod, "finalData.csv"))
allCensusedNegativeTrees <- analysesData[allCensusLiveTree=="yes" &
                                           positiveGrowthTree == "no" &
                                           BiomassGR > 0,]
nonAllCensusTrees <- analysesData[allCensusLiveTree == "no" & BiomassGR > 0]
origData <- analysesData[allCensusLiveTree == "yes" & positiveGrowthTree == "yes",]

for(indispecies in studySpecies){
  allHbestFormula <- allHbestFormulas[[indispecies]]
  theallHmodel <- allHbestModels[[indispecies]]
  speciesData <- origData[Species == indispecies,]
  speciesData1 <- allCensusedNegativeTrees[Species == indispecies,]
  speciesData1[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(speciesData$IniDBH)), 
                    Yearctd = Year-mean(speciesData$Year),
                    logHctd = log(H)-mean(log(speciesData$H)),
                    logSActd = log(IniFA+2.5)-mean(log(speciesData$IniFA+2.5)))]
  speciesData1$fittedlogY <- predict(theallHmodel, newdata = speciesData1,
                                   level = 0, se.fit = FALSE)
  R2_1 <- 1-sum((speciesData1$logY-speciesData1$fittedlogY)^2)/sum((speciesData1$logY-mean(speciesData1$logY))^2)
  speciesData2 <- nonAllCensusTrees[Species == indispecies,]
  speciesData2[,':='(logY = log(BiomassGR), 
                     logDBHctd = log(IniDBH)-mean(log(speciesData$IniDBH)), 
                     Yearctd = Year-mean(speciesData$Year),
                     logHctd = log(H)-mean(log(speciesData$H)),
                     logSActd = log(IniFA+2.5)-mean(log(speciesData$IniFA+2.5)))]
  speciesData2$fittedlogY <- predict(theallHmodel, newdata = speciesData2,
                                     level = 0, se.fit = FALSE)
  R2_2 <- 1-sum((speciesData2$logY-speciesData2$fittedlogY)^2)/sum((speciesData2$logY-mean(speciesData2$logY))^2)
  
  
  indioutput <- data.table(Species=indispecies, nonAllCensusR2 = R2_2, allCensusNegativeR2 = R2_1)
  if(indispecies == "All species"){
    alloutput <- indioutput
  } else {
    alloutput <- rbind(alloutput, indioutput)
  }
}

