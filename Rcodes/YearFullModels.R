rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table);library(MuMIn)
library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
selectionMethod <- "Year10Analyses"
analysesDataOrg <- fread(file.path(workPath, "data", selectionMethod, "finalData10.csv"))
analysesDataAll <- analysesDataOrg[allCensusLiveTree == "yes",]

studySpecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Minor species")



for(indispecies in studySpecies){
  speciesDataAll <- analysesDataAll[Species_Group == indispecies,]
  minABGR <- round(abs(min(speciesDataAll$BiomassGR)), 3)+0.01
  speciesDataAll[,':='(logY = log(BiomassGR+minABGR), 
                       logDBHctd = log(MidDBH)-mean(log(MidDBH)), 
                       Yearctd = MidYear-mean(MidYear),
                       logHctd = log(H)-mean(log(H)),
                       logSActd = log(MidFA)-mean(log(MidFA)))]
  if(indispecies == "All species"){
    fullModelAll <- lme(logY~logDBHctd+Yearctd+logHctd+logSActd+
                          logDBHctd:Yearctd+logDBHctd:logHctd+logDBHctd:logSActd+
                          Yearctd:logHctd+Yearctd:logSActd+
                          logHctd:logSActd,
                        data = speciesDataAll,
                        random = ~Yearctd+1|PlotID/uniTreeID, 
                        control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  } else {
    fullModelAll <- lme(logY~logDBHctd+Yearctd+logHctd+logSActd+
                          logDBHctd:Yearctd+logDBHctd:logHctd+logDBHctd:logSActd+
                          Yearctd:logHctd+Yearctd:logSActd+
                          logHctd:logSActd,
                        data = speciesDataAll,
                        random = ~Yearctd+1|PlotID/uniTreeID, 
                        control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
  }
  if(indispecies == "All species"){
    allfullModelsAll <- list(fullModelAll)
    names(allfullModelsAll) <- "All species"
  } else {
    allfullModelsAll[[indispecies]] <- fullModelAll
  }
  cat("Species", indispecies, "is done. \n")
  rm(fullModelAll)
}

fixedCoeffAll <- lapply(allfullModelsAll, function(x){
  data.table(summary(x)$tTable, keep.rownames = TRUE)[, ':='(marR2 = r.squaredGLMM(x)[1],
                                                             conR2 = r.squaredGLMM(x)[2])]})

for(i in 1:length(fixedCoeffAll)){
  if(i == 1){
    allcoeffs <- fixedCoeffAll[[1]][, Species:=names(fixedCoeffAll)[i]]
  } else {
    allcoeffs <- rbind(allcoeffs, fixedCoeffAll[[i]][, Species:=names(fixedCoeffAll)[i]])
  }
}

rm(indispecies,  i, minABGR,
   speciesDataAll)
save.image(file.path(workPath, "data", selectionMethod,
                     "FullYearModels.RData"))




