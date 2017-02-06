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
for(indispecies in studySpecies){
  allHbestFormula <- allHbestFormulas[[indispecies]]
  theallHmodel <- allHbestModels[[indispecies]]
  speciesData <- analysesData[Species == indispecies,]
  speciesData[,':='(logY = log(BiomassGR), 
                    logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                    Yearctd = 0,
                    logHctd = log(H)-mean(log(H)),
                    logSActd = log(IniFA+2.5)-mean(log(IniFA+2.5)))]
  speciesData$fittelogY <- predict(theallHmodel, newdata = speciesData,
                                   level = 0, se.fit = FALSE)
  plotRandom <- as.data.table(random.effects(theallHmodel)$PlotID,
                              keep.rownames = TRUE)
  names(plotRandom) <- c("PlotID", "PlotEffect")
  treeRandom <- as.data.table(random.effects(theallHmodel)$uniTreeID,
                              keep.rownames = TRUE)
  names(treeRandom) <- c("uniTreeID", "TreeEffect")
  treeRandom[, uniTreeID:=unlist(lapply(uniTreeID, 
                                        function(s) unlist(strsplit(s, split = "/", fixed = T))[2]))]
  speciesData <- setkey(speciesData, PlotID)[setkey(plotRandom, PlotID),
                                             nomatch = 0]
  speciesData <- setkey(speciesData, uniTreeID)[setkey(treeRandom, uniTreeID),
                                                nomatch = 0]
  
  speciesData[, logYResiduals:=logY-fittelogY-PlotEffect-TreeEffect]
  
  selectedplots <- unique(speciesData$PlotID)
  for(indiplot in selectedplots){
    indiplotdata <- speciesData[PlotID == indiplot,]
    linearRegression <- lm(logYResiduals~Year+logSActd:Year+Year:logHctd+Year:logDBHctd, data = indiplotdata)
    plotcoeffs <- as.data.table(summary(linearRegression)$coefficients,
                                keep.rownames = TRUE)
    tempplotcoeffs <- plotcoeffs[rn == "Year",.(Species = indispecies,
                                                NofTree = length(unique(indiplotdata$uniTreeID)),
                                                minYear = min(indiplotdata$IniYear),
                                                minSA = min(indiplotdata$IniFA),
                                                MeasurementLength = max(indiplotdata$IniFA)-
                                                  min(indiplotdata$IniFA),
                                            PlotID = indiplot, 
                                            YearEffect = Estimate,
                                            P = `Pr(>|t|)`)]
    if(nrow(plotcoeffs[rn == "Year:logHctd",])>0){
      tempplotcoeffs[, Interaction := plotcoeffs[rn == "Year:logHctd",]$Estimate]
    } else {
      tempplotcoeffs[, Interaction := 0]
    }
    if(indispecies == studySpecies[1] & indiplot == selectedplots[1]){
      allslopes <- tempplotcoeffs
    } else {
      allslopes <- rbind(allslopes, tempplotcoeffs)
    }
  }
}


for(indispecies in studySpecies){
  minYearModel_Year <- as.data.table(summary(lm(YearEffect~minYear, 
                                                data = allslopes[Species == indispecies]))$coefficients,
                                     keep.rownames = TRUE)
  MeasurementLModel_Year <- as.data.table(summary(lm(YearEffect~MeasurementLength, 
                                              data = allslopes[Species == indispecies]))$coefficients,
                                   keep.rownames = TRUE)
  
  YearEffect <- rbind(minYearModel_Year, MeasurementLModel_Year)
  YearEffect <- YearEffect[rn %in% c("NofTree", "minYear", "minSA", "MeasurementLength"),][,.(Variable = rn,
                                                                         YearCorrelation = paste(round(Estimate, 3),
                                                                                             "(P=",round(`Pr(>|t|)`, 4),")",
                                                                                             sep = ""))]
  minYearModel_Interaction <- as.data.table(summary(lm(Interaction~minYear, 
                                                data = allslopes[Species == indispecies]))$coefficients,
                                     keep.rownames = TRUE)
  MeasurementLModel_Interaction <- as.data.table(summary(lm(Interaction~MeasurementLength, 
                                                     data = allslopes[Species == indispecies]))$coefficients,
                                          keep.rownames = TRUE)
  
  Interaction <- rbind(minYearModel_Interaction,
                       MeasurementLModel_Interaction)
  Interaction <- Interaction[rn %in% c("NofTree", "minYear", "minSA", "MeasurementLength"),][,.(Variable = rn,
                                                                           InteractionCorrelation = paste(round(Estimate*(10^6), 2),
                                                                                                          "*10^-6 (P=",round(`Pr(>|t|)`, 4),")",
                                                                                                          sep = ""))]
  alleffects <- setkey(YearEffect, Variable)[setkey(Interaction, Variable),
                                             nomatch = 0]
  alleffects[,Species:=indispecies]
  if(indispecies == "All species"){
    alloutput <- alleffects
  } else {
    alloutput <- rbind(alloutput, alleffects)
  }

}

write.csv(alloutput, file.path(workPath, "TablesFigures",
                               "SampleStrategyEffects.csv"),
          row.names = FALSE)


