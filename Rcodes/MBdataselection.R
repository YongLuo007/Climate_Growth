rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- "J:/MBgrowth"
}

load(file.path(workPath, "data", "MBPSP.RData"))
names(allPSP) <- c("PlotID", "FMU", "TWP", "RGE", "PlotNumber",
                   "StandType", "StructureType", "VegeType",
                   "SStructure", "Moist", "YearEstabish", "Easting",
                   "Northing", "TreeNumber", "Species", "Distance",
                   "Angle", "DBH", "Height", "Status", "Class",
                   "TreeAge", "Year", "PlotSize")
orgData <- data.frame(allPSP) %>% data.table
print(length(unique(allPSP$PlotID))) # 425
# plot level selection
# 1. Natureal regeneration
allPSP <- allPSP[StandType == "NAT",]
print(length(unique(allPSP$PlotID))) # 271
# 2. minimum 3 times
allPSP[,measureTime:=length(unique(Year)), by = PlotID]
allPSP <- allPSP[measureTime>=4,]
print(length(unique(allPSP$PlotID))) # 172
# 3. plot size bigger than 500m^2
allPSP[,plotsizelength := length(unique(PlotSize)), by = PlotID]
allPSP[plotsizelength == 2,.(maxDistance=max(Distance)), by = PlotID] 
# PlotID maxDistance
# 1: 61-295       12.55
#*unique(allPSP[plotsizelength == 2,]$PlotSize)
# [1] 500.34 500.31
allPSP[plotsizelength == 2, PlotSize:=500.34]
#*unique(allPSP$PlotSize)
# [1] 500.34  50.01 250.17
allPSP <- allPSP[PlotSize>=500,]
print(length(unique(allPSP$PlotID))) # 169
set(allPSP, , c("measureTime", "plotsizelength"), NULL)

# 4 locations are available
range(allPSP$Easting) #[1] 325773 769020
range(allPSP$Northing) #[1] 5437011 6318185

# 5. FA information
allPSP[, ':='(baseTreeAge= TreeAge-Year+min(Year), baseYear=min(Year)), by = PlotID]
allPSP[,baseFA:=round(mean(baseTreeAge, na.rm = TRUE)), by = PlotID]
allPSP[,FA:=Year-baseYear+baseFA]
#* range(allPSP$FA)
# [1]   5 168
allPSP[,':='(baseTreeAge = NULL, baseYear = NULL, baseFA = NULL)]
allPSP[,':='(uniTreeID = paste(PlotID, "_", TreeNumber, sep = ""))]

# 6. plot level minimum 3 5-periods
selectedPlots <- unique(allPSP[,.(PlotID, Year)], by = c("PlotID", "Year"))

selectedPlots <- selectedPlots[order(PlotID, Year),]

newcol1 <- "MidYear"
newcol2 <- "FinYear"
col <- "Year"
selectedPlots[, c(newcol1, newcol2):=c(shift(.SD, n = 1, fill = NA, type = "lead"),
                                       shift(.SD, n = 2, fill = NA, type = "lead")),
              .SDcols = col, by = PlotID]
selectedPlots <- selectedPlots[!is.na(FinYear),]
selectedPlots[, ':='(first5P = MidYear - Year,
                     second5P = FinYear - MidYear)]
selectedPlots <- selectedPlots[first5P == 5 & second5P == 5,]
length(unique(selectedPlots$PlotID)) # 125 plots


# save a file for calculating compeitition
write.csv(allPSP, file.path(workPath, "data", "Year10Analyses", "forcompetitionIndex.csv"), row.names = FALSE)


allPSP <- allPSP[PlotID %in% unique(selectedPlots$PlotID),]
####$$$$%%%%%%%%%% tree level selection
# 1. select all alive trees
length(unique(allPSP$uniTreeID)) # 40554 trees
unique(allPSP$Status) #  1  2  4  0  9 NA  8
allPSP[,':='(unhealthyTrees=max(Status)), by = uniTreeID]
length(unique(allPSP[is.na(unhealthyTrees),]$uniTreeID)) # 20 trees that have NA status
allPSP <- allPSP[!is.na(unhealthyTrees),]
length(unique(allPSP$uniTreeID)) # 40534
length(unique(allPSP[unhealthyTrees==2,]$uniTreeID)) # 3233 trees have physical damage
allPSP <- allPSP[unhealthyTrees!=2,]
length(unique(allPSP$uniTreeID)) # 37301
length(unique(allPSP[unhealthyTrees==4,]$uniTreeID)) # 10 trees have severe insect attack
allPSP <- allPSP[unhealthyTrees!=4,]
length(unique(allPSP$uniTreeID)) # 37291
length(unique(allPSP[unhealthyTrees==9,]$uniTreeID)) # 306 trees have unknown causes of death
allPSP <- allPSP[unhealthyTrees!=9,]
length(unique(allPSP$uniTreeID)) # 36985
length(unique(allPSP[unhealthyTrees==8,]$uniTreeID)) # 1 tree's death due to wind/snow 
allPSP <- allPSP[unhealthyTrees!=8,]
length(unique(allPSP$uniTreeID)) # 36984

###############################
# there are 36984 trees
###############################
allPSP[,':='(firstPlotYear = min(Year), lastPlotYear = max(Year)), by = PlotID]
allPSP[,':='(firstTreeYear = min(Year), lastTreeYear = max(Year)), by = uniTreeID]
allPSP[firstPlotYear == firstTreeYear, RecruitMent := 0]
allPSP[firstPlotYear != firstTreeYear, RecruitMent := 1]
length(unique(allPSP[RecruitMent == 1,]$uniTreeID)) # 9 trees recruited (interesting)
allPSP[lastPlotYear == lastTreeYear, Mortality := 0]
allPSP[lastPlotYear != lastTreeYear, Mortality := 1]
length(unique(allPSP[Mortality == 1,]$uniTreeID)) # 16297 trees (make sense)



allPSP <- allPSP[order(uniTreeID, Year),]
shiftcol <- c("Year", "DBH", "Height", 
              "uniTreeID", "FA")
midcol <- c("MidYear", "MidDBH", "MidHeight",  "MiduniTreeID",
            "MidFA") 

fincol <- c("FinYear", "FinDBH", "FinHeight",  "tempuniTreeID",
                     "FinFA") 

allPSP[, (midcol) := data.table::shift(x = .SD,
                                       n = 1, fill = NA, type = "lead"), .SDcols = shiftcol,
       by = uniTreeID]
allPSP[, (fincol) := data.table::shift(x = .SD,
                                       n = 2, fill = NA, type = "lead"), .SDcols = shiftcol,
       by = uniTreeID]


simplePSP <- allPSP[!is.na(tempuniTreeID), ][,c("MiduniTreeID", "tempuniTreeID") := NULL]

length(unique(simplePSP$uniTreeID)) # 30295
simplePSP <- simplePSP[,.(PlotID, uniTreeID, Species,IniYear = Year, IniFA = FA,
                       IniDBH = DBH, IniHeight = Height, 
                       MidYear, MidDBH, MidHeight, MidFA,
                       FinYear, FinDBH, FinHeight, FinFA)]



inputData <- data.table(data.frame(simplePSP))
# "AS" 
inputData[Species == "AS",species:="black spruce"]
# "BA"
inputData[Species == "BA",species:="balsam poplar"]
# "BF" 
inputData[Species == "BF",species:="balsam fir"]
# "BO"
inputData[Species == "BO",species:="white oak"] # bur oak to
# "BS"
inputData[Species == "BS",species:="black spruce"]
# "EC"
inputData[Species == "EC",species:="western redcedar"] # cedar
# "JP"
inputData[Species == "JP",species:="jack pine"]
# "MM"
inputData[Species == "MM",species:="silver maple"]
# "RP"
inputData[Species == "RP",species:="red pine"]
# "TA"
inputData[Species == "TA",species:="trembling aspen"]
# "TL"
inputData[Species == "TL",species:="tamarack larch"]
# "WB"
inputData[Species == "WB",species:="white birch"]
# "WE"
inputData[Species == "WE",species:="white elm"]
# "WS"
inputData[Species == "WS",species:="white spruce"]

inputData[,length:=FinYear-IniYear]

inputData <- inputData[length == 10,]
inputData[,measureLength:=length(unique(IniYear)), by = PlotID]
inputData <- inputData[measureLength>=2,]
unique(inputData$length) # 10
length(unique(inputData$PlotID)) # 125 plots
# recruitment trees
inputData[, ':='(plotfirstY = min(IniYear), plotlastY = max(FinYear)), by = PlotID]
inputData[,':='(treefirstY = min(IniYear), treelastY = max(FinYear)), by = uniTreeID]
inputData[,':='(regTree = 0, deadTree = 0)]
inputData[plotfirstY != treefirstY, regTree:=1]
inputData[plotlastY != treelastY, deadTree:=1]
inputData[,allCensusLiveTree := "no"]
inputData[regTree == 0 & deadTree == 0, allCensusLiveTree := "yes"]
inputData[,positiveGrowth:="yes"]
inputData[FinDBH<=IniDBH, positiveGrowth:="no"]
inputData[,lengthNegativeG:=length(unique(positiveGrowth)), by = uniTreeID]
inputData[,positiveGrowthtemp := "no"]
inputData[positiveGrowth == "yes" & lengthNegativeG == 1, positiveGrowthtemp := "yes"]

inputData <- inputData[,.(PlotID, uniTreeID, allCensusLiveTree, 
                          positiveGrowthTree = positiveGrowthtemp, 
                          Species_Org = Species, Species_Std = species,
                          IniYear, IniFA, IniDBH, MidYear, MidFA, MidDBH, 
                          FinYear, FinDBH, FinFA)]
selectedplotsummary <- unique(inputData[,.(PlotID, IniYear, FinYear)], by = c("PlotID", "IniYear"))
write.csv(selectedplotsummary, file.path(workPath, "data", "Year10Analyses", "plotsummary.csv"),row.names=FALSE)


selectedPSP <- allPSP[PlotID %in% unique(inputData$PlotID), ]
selectedPSPMasterTable <- unique(selectedPSP[,.(PlotID, FMU, TWP, RGE, StandType, StructureType,
                                                VegeType, VegeType, Easting, Northing)],
                                 by = "PlotID")
write.csv(selectedPSPMasterTable, file.path(workPath, "data", "Year10Analyses", "selectedPlotMasterTable.csv"),row.names=FALSE)

write.csv(inputData, file.path(workPath, "data", "Year10Analyses",
                               "MBdataSimplified.csv"),
          row.names = FALSE)




