rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table);library(MuMIn)
library(parallel)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "newAllDataRescaledComp.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
analysesData <- analysesData[allCensusLiveTree == "yes" & positiveGrowthTree == "yes",]

unique(analysesData$species)
# [1] "jack pine"        "trembling aspen"  "black spruce"    
# [4] "balsam poplar"    "tamarack larch"   "white spruce"    
# [7] "white birch"      "balsam fir"       "western redcedar"
# [10] "white oak"        "red pine"         "white elm"   

#### assign wood density by species
# source
# http://cfs.nrcan.gc.ca/pubwarehouse/pdfs/11744.pdf
# ovendry wood density is used
# wood density unit: kg/m3
analysesData[species == "balsam fir", woodDensity := 367]
analysesData[species == "tamarack larch", woodDensity := 530]
analysesData[species == "jack pine", woodDensity := 451]
analysesData[species == "red pine", woodDensity := 419]
analysesData[species == "black spruce", woodDensity := 457]
analysesData[species == "white spruce", woodDensity := 404]
analysesData[species == "white birch", woodDensity := 607]
analysesData[species == "white elm", woodDensity := 617]
analysesData[species == "white oak", woodDensity := 775]
analysesData[species == "balsam poplar", woodDensity := 409]
analysesData[species == "trembling aspen", woodDensity := 424]
analysesData[species == "western redcedar", woodDensity := 338]


a <- analysesData[,.(meanBiomassGR = mean(BiomassGR), meanBAGR = mean(BAGR),
                     meanWD = mean(woodDensity)), by = species]

analysesData[,':='(logY = log(BiomassGR), 
                  logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                  Yearctd = Year-mean(Year),
                  logIntraHctd = log(IntraH+1)-mean(log(IntraH+1)),
                  logInterHctd = log(InterH+1)-mean(log(InterH+1)),
                  logSActd = log(IniFA+2.5)-mean(log(IniFA+2.5)),
                  logWD = log(woodDensity)-mean(log(woodDensity)))]
source(file.path(workPath, "Rcodes", "Rfunctions", "mixedModelSelection.R"))
tempoutput <- mixedModelSelection(DV = "logY", 
                                  IDV = c("logDBHctd", "Yearctd", "logIntraHctd", 
                                          "logInterHctd", "logSActd"),
                                  maxInteraction = 2,
                                  ICTerm = "AIC",
                                  ICCut = 2,
                                  data = analysesData,
                                  random = ~1|uniTreeID, 
                                  control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
bestFormu <- as.formula(paste("logY~", paste(tempoutput$bestIDV, collapse = "+")))
bestModel <- lme(fixed = bestFormu,
                 data = analysesData,
                 random = ~1|uniTreeID, 
                 control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))

