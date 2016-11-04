rm(list=ls())
library(data.table)
studySpecies <- c("All species", "Jack pine", 
                  "Trembling aspen", "Black spruce", "Other species")
b5 <- list()
a <- list()
for(indispecies in studySpecies){
  if(indispecies == "All species"){
    addon <- "_2"
  } else if(indispecies == "Jack pine"){
    addon <- "_6"
  } else if(indispecies == "Trembling aspen"){
    addon <- "_6"
  } else if(indispecies == "Black spruce") {
    addon <- "_6"
  } else if(indispecies == "Other species"){
    addon <- "_6"
  }
  for(k in 1:3){
    workPath <- "~/GitHub/Climate_Growth"
    load(file.path(workPath, "data", "MCMC analyses", 
                   indispecies, 
                   paste("chain", k, addon,
                         ".RData", sep = "")))
    if(k == 1){
      allcoda <- get(paste("coda", k, sep = ""))
    } else {
      allcoda[k] <- get(paste("coda", k, sep = ""))
    }
  }
  rm(i, coda1, coda2, coda3, k)
  codatable <- data.table(as.matrix(allcoda))[,':='(time = time(allcoda), 
                                            chain = sort(rep(1:3, 2500)))]
  b5[[indispecies]] <- ggplot(data = codatable, aes(x = time, y = b5))+
    geom_line(aes(group = chain, col = as.factor(chain)))
  a[[indispecies]] <- ggplot(data = codatable, aes(x = time, y = a))+
    geom_line(aes(group = chain, col = as.factor(chain)))
   # check the convergency
  gelmantest <- gelman.diag(allcoda, multivariate = FALSE)
  gelmanResults <- data.table(gelmantest$psrf, keep.rownames = TRUE)
  # the main coefficients
  cat("For", indispecies, "\n")
  print(gelmanResults[rn %in% c("a", paste("b", 1:5, sep = ""), 
                                paste("c", 1:10, sep = ""),
                                paste("d", 1:2, sep = "")), ])
 # summary the posterior distributions
  allPosterialDisSumm <- summary(allcoda)
  allPosterialDisSumm1 <- data.table(cbind(allPosterialDisSumm$statistics, allPosterialDisSumm$quantiles),
                                    keep.rownames = TRUE)
  fixedEffect_indispecies <- allPosterialDisSumm1[rn %in% c("a", paste("b", 1:5, sep = ""), paste("c", 1:10, sep = ""),
                                                paste("d", 1:2, sep = "")), ][
                                                  ,.(Species = indispecies, Variable = rn, Mean = Mean, 
                                                     SD = SD, Lower95 = `2.5%`, Upper = `97.5%`)]
  OverAllTrend <- allPosterialDisSumm1[rn %in% c(paste("predictY1[", 1:50, "]", sep = "")),]
  TrendWithRBI <- allPosterialDisSumm1[rn %in% c(paste("predictY2[", 1:500, "]", sep = "")),]
  Figure2Data_indispecies <- rbind(OverAllTrend[,.(Species = indispecies, Direction = "Overall Trend", 
                                                   Year = YearPredict+mean(speciesData$Year), RBI = 0,
                                                   PredictedABGR = Mean, PredictedABGR_Lower = `2.5%`,
                                                   PredictedABGR_Upper = `97.5%`)],
                                   TrendWithRBI[,.(Species = indispecies, Direction = "Change with RBI",
                                                   Year = YearRBITable$Year, RBI = YearRBITable$RBI,
                                                   PredictedABGR = Mean, PredictedABGR_Lower = `2.5%`,
                                                   PredictedABGR_Upper = `97.5%`)])
  IntraHTrendwithRBI <- allPosterialDisSumm1[rn %in% paste("PredictIntraH[", 1:500, "]", sep = ""),]
  InterHTrendwithRBI <- allPosterialDisSumm1[rn %in% paste("PredictInterH[", 1:500, "]", sep = ""),]
  Figure3Data_indispecies <- rbind(IntraHTrendwithRBI[,.(Direction = "IntraH", Competition = "IntraH",
                                                         Species = indispecies, Year = YearRBITable$Year,
                                                         RBI = YearRBITable$RBI, Value = Mean, 
                                                         Value_Lower = `2.5%`, Value_Upper = `97.5%`,
                                                         Main = 0)],
                                   InterHTrendwithRBI[,.(Direction = "InterH", Competition = "InterH",
                                                         Species = indispecies, Year = YearRBITable$Year,
                                                         RBI = YearRBITable$RBI, Value = Mean, 
                                                         Value_Lower = `2.5%`, Value_Upper = `97.5%`,
                                                         Main = 0)],
                                   data.table(Direction = c("IntraH", "InterH"), 
                                              Competition = c("IntraH", "InterH"),
                                              Species = indispecies, Year = mean(speciesData$Year),
                                              RBI = 0, 
                                              Value = allPosterialDisSumm1[rn %in% c("b2", "b3"), ]$Mean,
                                              Value_Lower = allPosterialDisSumm1[rn %in% c("b2", "b3"), ]$`2.5%`,
                                              Value_Upper = allPosterialDisSumm1[rn %in% c("b2", "b3"), ]$`97.5%`,
                                              Main = 1))
  # summary the posterior distribution for main Year effect (b5) and its change with RBI (c10)
  allcodaDataTable <- data.table(as.matrix(allcoda))[,.(b2, b3, b5, c10, d1, d2)]
  
  tempdensity <- density(allcodaDataTable$b5, n = 500,
                        from = min(allcodaDataTable$b5),
                        to = max(allcodaDataTable$b5))
  Yeareffect_indispecies <- data.table(Species = indispecies,
                                       Yeareffect = tempdensity$x,
                                       Density = tempdensity$y)
  rm(tempdensity)
  tempdensity <- density(allcodaDataTable$c10, n = 500,
                         from = min(allcodaDataTable$c10),
                         to = max(allcodaDataTable$c10))
  YeareffectWithRBI_indispecies <- data.table(Species = indispecies,
                                       Yeareffect = tempdensity$x,
                                       Density = tempdensity$y)
  if(indispecies == "All species"){
    fixedEffect <- fixedEffect_indispecies
    Figure2Data <- Figure2Data_indispecies
    Figure3Data <- Figure3Data_indispecies
    Yeareffect <- Yeareffect_indispecies
    YeareffectWithRBI <- YeareffectWithRBI_indispecies
  } else {
    fixedEffect <- rbind(fixedEffect, fixedEffect_indispecies)
    Figure2Data <- rbind(Figure2Data, Figure2Data_indispecies)
    Figure3Data <- rbind(Figure3Data, Figure3Data_indispecies)
    Yeareffect <- rbind(Yeareffect, Yeareffect_indispecies)
    YeareffectWithRBI <- rbind(YeareffectWithRBI, YeareffectWithRBI_indispecies)
  }
}
rm(coda1, coda2, coda3, allPosterialDisSumm, allPosterialDisSumm1,
   chain1, chain2,chain3, data, inits, OverAllTrend, TrendWithRBI, speciesData,
   tempdensity, workPath, Yeareffect_indispecies, YeareffectWithRBI_indispecies,
   Figure2Data_indispecies, Figure3Data_indispecies, fixedEffect_indispecies)
workPath <- "~/GitHub/Climate_Growth"
save.image(file.path(workPath, "data", "Bayesian_YearModels.RData"))
