rm(list=ls())
library(data.table)
studySpecies <- c("All", "JP", "TA", "BS", "Other")
b2 <- list()
a <- list()
library(coda)
workPath <- "~/GitHub/Climate_Growth"
for(indispecies in studySpecies[2:5]){
  if(indispecies == "All"){
    addon <- ""
  } else if(indispecies == "JP"){
    addon <- "_3"
  } else if(indispecies == "TA"){
    addon <- "_3"
  } else if(indispecies == "BS") {
    addon <- "_2"
  } else if(indispecies == "Other"){
    addon <- "_3"
  }
  for(k in 1:3){
    load(file.path(workPath, "Results", "JAGS", "YearReduced",
                   indispecies, 
                   paste("chain", k, addon,
                         ".RData", sep = "")))
    if(k == 1){
      allcoda <- get(paste("coda", k, sep = ""))
    } else {
      allcoda[k] <- get(paste("coda", k, sep = ""))
    }
  }
  rm(coda1, coda2, coda3, k)
  codatable <- data.table(as.matrix(allcoda))[,':='(time = time(allcoda), 
                                            chain = sort(rep(1:3, length(time(allcoda)))))]
  colNum <- which(names(codatable) %in% c("a", paste("b", 1:7, sep = ""), 
                                           paste("c", 1:21, sep = "")))
  subcoda <- allcoda[,colNum, drop = TRUE]
   # check the convergency
  gelmantest <- gelman.diag(subcoda, multivariate = FALSE)
  gelmanResults <- data.table(gelmantest$psrf, keep.rownames = TRUE)
  # the main coefficients
  cat("For", indispecies, "\n")
  print(gelmanResults[rn %in% c("a", paste("b", 1:7, sep = ""), 
                                paste("c", 1:21, sep = "")),])
 # summary the posterior distributions
  allPosterialDisSumm <- summary(allcoda)
  allPosterialDisSumm1 <- data.table(cbind(allPosterialDisSumm$statistics, allPosterialDisSumm$quantiles),
                                    keep.rownames = TRUE)
  fixedEffect_indispecies <- allPosterialDisSumm1[,.(Species = indispecies, rn, Mean = Mean, 
                                                     SD = SD, Lower95 = `2.5%`, Upper95 = `97.5%`)]
  fixedEffect_indispecies <- setkey(fixedEffect_indispecies, rn)[setkey(speciesinits[,.(rn = coeffs, Variable)], rn),
                                                                 nomatch = 0][,.(Species, coeffs = rn, Variable,
                                                                                 Mean, SD, Lower95, Upper95)]
  
  OverAllTrend <- allPosterialDisSumm1[rn %in% c(paste("predictBAGR1[", 1:50, "]", sep = "")),]
  if(nrow(OverAllTrend)>0){
    OverAllTrend <- OverAllTrend[,.(Species = indispecies, Direction = "Overall Trend", 
                                                   Year = YearPredict+mean(speciesData$Year), RBI = 0,
                                                   PredictedABGR = Mean, PredictedABGR_Lower = `2.5%`,
                                                   PredictedABGR_Upper = `97.5%`)]
  } else {
    OverAllTrend <- data.table(Species = character(), Direction = character(),
                               Year = numeric(), RBI = numeric(), PredictedABGR = numeric(),
                               PredictedABGR_Lower = numeric(), PredictedABGR_Upper = numeric())
  }
  TrendWithRBI <- allPosterialDisSumm1[rn %in% c(paste("predictBAGR2[", 1:500, "]", sep = "")),]
  if(nrow(TrendWithRBI)>0){
    TrendWithRBI <- TrendWithRBI[,.(Species = indispecies, Direction = "Change with RBI",
                                                   Year = YearRBITable$Year, RBI = YearRBITable$RBI,
                                                   PredictedABGR = Mean, PredictedABGR_Lower = `2.5%`,
                                                  PredictedABGR_Upper = `97.5%`)]
  } else {
    TrendWithRBI <- data.table(Species = character(), Direction = character(),
                               Year = numeric(), RBI = numeric(), PredictedABGR = numeric(),
                               PredictedABGR_Lower = numeric(), PredictedABGR_Upper = numeric())
  }
   Figure2Data_indispecies <- rbind(OverAllTrend, TrendWithRBI) 
 
  
  IntraHTrendwithRBI <- allPosterialDisSumm1[rn %in% paste("predictIntraHeffect[", 1:500, "]", sep = ""),]
  if(nrow(IntraHTrendwithRBI)>0){
    IntraHTrendwithRBI <- IntraHTrendwithRBI[,.(Direction = "IntraH", Competition = "IntraH",
                                                Species = indispecies, 
                                                Year = sort(unique(YearRBITable$Year)),
                                                RBI = 0, Value = Mean, 
                                                Value_Lower = `2.5%`, Value_Upper = `97.5%`,
                                                Main = 0)]
  } else {
    IntraHTrendwithRBI <- data.table(Direction = character(), Competition = character(),
                                     Species = character(), Year = numeric(), RBI = numeric(),
                                     Value = numeric(), Value_Lower = numeric(), Value_Upper = numeric(),
                                     Main = numeric())
  }
  InterHTrendwithRBI <- allPosterialDisSumm1[rn %in% paste("predictInterHeffect[", 1:500, "]", sep = ""),]
  if(nrow(InterHTrendwithRBI)>0){
    InterHTrendwithRBI <- InterHTrendwithRBI[,.(Direction = "InterH", Competition = "InterH",
                                                Species = indispecies, 
                                                Year = sort(unique(YearRBITable$Year)),
                                                RBI = 0, Value = Mean, 
                                                Value_Lower = `2.5%`, Value_Upper = `97.5%`,
                                                Main = 0)]
  } else {
    InterHTrendwithRBI <- data.table(Direction = character(), Competition = character(),
                                     Species = character(), Year = numeric(), RBI = numeric(),
                                     Value = numeric(), Value_Lower = numeric(), Value_Upper = numeric(),
                                     Main = numeric())
  }
  mainCompetitioneff <- allPosterialDisSumm1[rn %in% c("b3", "b4"), ][
    , ':='(Species = indispecies, Year = mean(speciesData$Year), RBI = 0)]
  mainCompetitioneff[rn == "b3", Competition := "IntraH"]
  mainCompetitioneff[rn == "b4", Competition := "InterH"]
  mainCompetitioneff <- mainCompetitioneff[,.(Direction = Competition, Competition, Species, Year, RBI, 
                                              Value = Mean, Value_Lower = `2.5%`, Value_Upper = `97.5%`,
                                              Main = 1)]

  Figure3Data_indispecies <- rbind(IntraHTrendwithRBI, InterHTrendwithRBI, mainCompetitioneff)
  # # summary the posterior distribution for main Year effect (b2) and its change with RBI (c10)
  # allcodaDataTable <- data.table(as.matrix(allcoda))[,.(b2, b3, b2, c10, d1, d2)]
  # 
  # tempdensity <- density(allcodaDataTable$b2, n = 500,
  #                       from = min(allcodaDataTable$b2),
  #                       to = max(allcodaDataTable$b2))
  # Yeareffect_indispecies <- data.table(Species = indispecies,
  #                                      Yeareffect = tempdensity$x,
  #                                      Density = tempdensity$y)
  # rm(tempdensity)
  # tempdensity <- density(allcodaDataTable$c10, n = 500,
  #                        from = min(allcodaDataTable$c10),
  #                        to = max(allcodaDataTable$c10))
  # YeareffectWithRBI_indispecies <- data.table(Species = indispecies,
  #                                      Yeareffect = tempdensity$x,
  #                                      Density = tempdensity$y)
  if(indispecies == "JP"){
    fixedEffect <- fixedEffect_indispecies
    Figure2Data <- Figure2Data_indispecies
    Figure3Data <- Figure3Data_indispecies
    # Yeareffect <- Yeareffect_indispecies
    # YeareffectWithRBI <- YeareffectWithRBI_indispecies
  } else {
    fixedEffect <- rbind(fixedEffect, fixedEffect_indispecies)
     Figure2Data <- rbind(Figure2Data, Figure2Data_indispecies)
    Figure3Data <- rbind(Figure3Data, Figure3Data_indispecies)
    # Yeareffect <- rbind(Yeareffect, Yeareffect_indispecies)
    # YeareffectWithRBI <- rbind(YeareffectWithRBI, YeareffectWithRBI_indispecies)
  }
}

rm(coda1, coda2, coda3, allPosterialDisSumm, allPosterialDisSumm1,
   chain1, chain2,chain3, data, inits, OverAllTrend, TrendWithRBI, speciesData,
   tempdensity, workPath, Yeareffect_indispecies, YeareffectWithRBI_indispecies,
   Figure2Data_indispecies, Figure3Data_indispecies, fixedEffect_indispecies)
workPath <- "~/GitHub/Climate_Growth"
save.image(file.path(workPath, "data", "Bayesian_YearModels.RData"))
