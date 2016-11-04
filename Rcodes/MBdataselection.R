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
allPSP <- allPSP[measureTime>=3,]
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
# 4. FA information
allPSP[, ':='(baseTreeAge= TreeAge-Year+min(Year), baseYear=min(Year)), by = PlotID]
allPSP[,baseFA:=round(mean(baseTreeAge, na.rm = TRUE)), by = PlotID]
allPSP[,FA:=Year-baseYear+baseFA]
#* range(allPSP$FA)
# [1]   5 168
allPSP[,':='(baseTreeAge = NULL, baseYear = NULL, baseFA = NULL)]


allPSP[,':='(uniTreeID = paste(PlotID, "_", TreeNumber, sep = ""))]
firsttime <- FALSE
if(firsttime){
   #alive trees for competition calculation

   CIcompetitionData <- allPSP[Status == 1 | Status == 0,]
   CIcompetitionData[PlotID == "46-132", Distance:=Distance/100]
   CIcompetitionData <- CIcompetitionData[,.(PlotNumber = PlotID, TreeNumber = TreeNumber,
                                             Year, Distance, Angle, DBH, Species)]
   
   source(file.path(workPath, "Rcodes", "Rfunctions","HeghyiCICalculation.R"))
   CIdata <- HeghyiCICalculation(data = CIcompetitionData,
                                 maxRadius = 12.62,
                                 sizeIndex = "DBH")
   
   
   write.csv(CIdata[,.(PlotNumber, uniTreeID Year, Hegyi, IntraHegyiRatio)],
             ,row.names = F)
}
rm(firsttime)



####$$$$%%%%%%%%%% tree level selection
# 1. select all alive trees
length(unique(allPSP$uniTreeID)) # 57157 trees

unique(allPSP$Status) #  1  2  4  0  9 NA  8
allPSP[,':='(unhealthyTrees=max(Status)), by = uniTreeID]
length(unique(allPSP[is.na(unhealthyTrees),]$uniTreeID)) # 20 trees that have NA status
allPSP <- allPSP[!is.na(unhealthyTrees),]
length(unique(allPSP$uniTreeID)) # 57137
length(unique(allPSP[unhealthyTrees==2,]$uniTreeID)) # 4430 trees have physical damage
allPSP <- allPSP[unhealthyTrees!=2,]
length(unique(allPSP$uniTreeID)) # 52707
length(unique(allPSP[unhealthyTrees==4,]$uniTreeID)) # 10 trees have severe insect attack
allPSP <- allPSP[unhealthyTrees!=4,]
length(unique(allPSP$uniTreeID)) # 52697
length(unique(allPSP[unhealthyTrees==9,]$uniTreeID)) # 360 trees have unknown causes of death
allPSP <- allPSP[unhealthyTrees!=9,]
length(unique(allPSP$uniTreeID)) # 52337
length(unique(allPSP[unhealthyTrees==8,]$uniTreeID)) # 1 tree's death due to wind/snow 
allPSP <- allPSP[unhealthyTrees!=8,]
length(unique(allPSP$uniTreeID)) # 52336

###############################
# there are 52336 trees
###############################
allPSP[,':='(firstPlotYear = min(Year), lastPlotYear = max(Year)), by = PlotID]
allPSP[,':='(firstTreeYear = min(Year), lastTreeYear = max(Year)), by = uniTreeID]
allPSP[firstPlotYear == firstTreeYear, RecruitMent := 0]
allPSP[firstPlotYear != firstTreeYear, RecruitMent := 1]
length(unique(allPSP[RecruitMent == 1,]$uniTreeID)) # 11 trees recruited (interesting)
allPSP[lastPlotYear == lastTreeYear, Mortality := 0]
allPSP[lastPlotYear != lastTreeYear, Mortality := 1]
length(unique(allPSP[Mortality == 1,]$uniTreeID)) # 20627 trees (make sense)
# 20627/52336 = 0.3941 
# 20627/nrow(allPSP) = 0.11 ()

# 2. measured at least two times so that the growth rate can be derive from
allPSP[,':='(maxLength = length(DBH)),
       by = uniTreeID]
allPSP <- allPSP[maxLength >= 2,]
length(unique(allPSP$uniTreeID)) # 45862
###############################
# there are 45862 trees
###############################


allPSP <- allPSP[order(uniTreeID, Year),]
allPSP <- allPSP[,c("FinYear", "FinDBH", "FinHeight",  "tempuniTreeID",
                     "FinFA") 
                  := shift(x = allPSP[,.(Year, DBH, Height, uniTreeID,
                                         FA)],
                           n = 1, fill = NA, type = "lead")]
allPSP <- allPSP[tempuniTreeID == uniTreeID,][,tempuniTreeID := NULL]
allPSP <- allPSP[Year != FinYear,]
length(unique(allPSP$uniTreeID)) # 45861
simplePSP <- allPSP[,.(PlotID, uniTreeID, Species,IniYear = Year, IniFA = FA,
                       IniDBH = DBH, IniHeight = Height, FinYear, 
                       FinDBH, FinHeight, FinFA)]

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


write.csv(inputData, file.path(workPath, "data","MBdataSimplified.csv"),
          row.names = FALSE)




