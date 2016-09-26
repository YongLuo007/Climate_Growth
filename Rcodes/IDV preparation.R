rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
thedata <- read.csv(file.path(workPath, "data", "MBdataSimplified.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>%
  data.table

source("~/GitHub/landwebNRV/landwebNRV/R/biomassCalculation.R")
inputData$IniBiomass <- biomassCalculation(species = inputData$species,
                                           DBH = inputData$IniDBH, 
                                           height = inputData$IniHeight)$biomass
inputData$FinBiomass <- biomassCalculation(species = inputData$species,
                                           DBH = inputData$FinDBH,
                                           height = inputData$FinHeight)$biomass
inputData[,':='(IniBA = 3.1415926*((IniDBH/2)^2), 
                FinBA = 3.1415926*((FinDBH/2)^2))]
inputData[,':='(Year = (FinYear+IniYear)/2,
                BAGR = (FinBA-IniBA)/(FinYear-IniYear),
                BiomassGR = (FinBiomass - IniBiomass)/(FinYear-IniYear))]

thedata[,NofTree:=length(unique(uniTreeID)), by = c("Year", "PlotID")]
thedata[,MinNofTree:=min(NofTree), by = PlotID]
thedata <- thedata[MinNofTree>=20, ][,':='(NofTree = NULL, MinNofTree = NULL)]
# print(length(unique(thedata$PlotID))) # 148 plots

thedata[,':='(TreeMLength = length(IniDBH)),
             by = uniTreeID]
thedata[, ':='(PlotMLength = length(unique(IniYear))), by = PlotID]
thedata[,Status:="IngrowthOrDead"]
thedata[TreeMLength == PlotMLength, Status:="Survive"]
set(thedata, , c("TreeMLength", "PlotMLength"), NULL)
minDBHdata <- thedata[Status == "Survive",.(PlotID, IniYear, IniDBH, IniHeight)][
  ,.(minDBH=min(IniDBH)), by = c("PlotID", "IniYear")]
thedata <- setkey(thedata, PlotID, IniYear)[setkey(minDBHdata, PlotID, IniYear), 
                                                      nomatch = 0]
thedata[,':='(maxDBH=max(IniDBH)), by = c("PlotID", "IniYear")]
thedata[,Dominance_indiSize:=100*(IniDBH-minDBH)/(maxDBH-minDBH)]

minBiomassdata <- thedata[Status == "Survive",.(PlotID, IniYear, IniBiomass, IniHeight)][
  ,.(minBiomass=min(IniBiomass)), by = c("PlotID", "IniYear")]
thedata <- setkey(thedata, PlotID, IniYear)[setkey(minBiomassdata, PlotID, IniYear), 
                                            nomatch = 0]
thedata[,':='(maxBiomass=max(IniBiomass)), by = c("PlotID", "IniYear")]
thedata[,Dominance_indiBiomass:=100*(IniBiomass-minBiomass)/(maxBiomass-minBiomass)]

# Dominance_indiSize is the dominance index that is associated with each tree's size
# 0 to 100, 0 is the smallest tree and 100 is the biggest tree

source(file.path(workPath, "Rcodes", "Rfunctions", "dominanceAssigner.R"))
thedata <- dominanceAssigner(inputData = thedata, 
                                   dominanceCut = seq(0, 1, by = 0.1),
                                   dominanceMethod = "Biomass")
thedata <- dominanceAssigner(inputData = thedata, 
                             dominanceCut = seq(0, 1, by = 0.1),
                             dominanceMethod = "BA")

thedata[, minBAGR:=round(min(BAGR), 3), by = uniTreeID]
thedata <- thedata[Status == "Survive" & minBAGR > 0,]

climateinputs <- read.csv(file.path(workPath, "data", "plotClimates.csv"), header = TRUE,
                          stringsAsFactors = FALSE) %>%
  data.table
nrow(thedata) #61373
thedata <- setkey(thedata, PlotID, IniYear, FinYear)[setkey(climateinputs, PlotID, IniYear, FinYear),
                                                               nomatch = 0]

nrow(thedata) #61373




speciesTable <- unique(thedata[,.(Species, uniTreeID)], by = "uniTreeID")[
  , .(NumberofTree=length(uniTreeID)), by = c("Species")] 

speciesTable[,totNumbofTree:=sum(NumberofTree), by = Species]
speciesTable <- speciesTable[order(-totNumbofTree),.(Species,  NumberofTree)]

set(thedata, ,c("Status", "minDBH", "maxDBH", "minBiomass", "maxBiomass","minBAGR"), NULL)

write.csv(thedata,
          file.path(workPath, "data", "MBdatafinal.csv"),
          row.names = FALSE)


