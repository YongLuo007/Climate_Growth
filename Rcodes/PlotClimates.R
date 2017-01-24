rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- "J:/MBgrowth"
}

inputData <- read.csv(file.path(workPath, "data", "plotsummary.csv"), header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table
climateData <- read.csv(file.path(workPath, "data", "StudyAreaClimates_BiomSIM", "plotMonthlyClimates.csv"),
                        header = TRUE, stringsAsFactors = FALSE) %>% data.table

names(climateData) <- c("PlotID", "Year_Month", "Temp", "Prep", "PET")

climateData[,':='(Year = as.numeric(unlist(lapply(Year_Month, function(x) unlist(strsplit(x, "/", fixed = T))[1]))),
                  Month = as.numeric(unlist(lapply(Year_Month, function(x) unlist(strsplit(x, "/", fixed = T))[2]))))]
CO2data <- read.csv(file.path(workPath, "data", "StudyAreaClimates_BiomSIM", "CO2.csv"),
                        header = TRUE, stringsAsFactors = FALSE) %>% data.table

climateData <- setkey(climateData, Year, Month)[setkey(CO2data[,.(Year, Month, CO2)], Year, Month),
                                                nomatch = 0]

# crt is current that for the climates for the current growth year
# ie. measurement year 2005-2010, climate window 2005-2010, 5years
# pre is previous that for the climate that contains the previous year climate 
# ie., measurement year 2005-2010, climate window 2004-2010, 6years
# GS is growing season from May to Septempber, 
# while NONGS is not growing season from previous October to current April
# long term is defined as 1985 to 2011
climateData[,CMI:=Prep - PET]
climateData[Month>=10, Year:=Year+1]
AnnualClimate <- climateData[Year>=1985 & Year<= 2011,][,.(AT = mean(Temp), AP = sum(Prep),
                                                           APET = sum(PET),
                                                           ACMI = sum(CMI), CO2 = mean(CO2)), 
                                                        by = c("PlotID", "Year")]
 GSClimate <- climateData[Year>=1985 & Year<= 2011 & Month<=9 & Month >= 5,][
   ,.(GST = mean(Temp), GSP = sum(Prep), GSPET = sum(PET), GSCMI = sum(CMI), GSCO2 = mean(CO2)), by = c("PlotID", "Year")]
 
 NONGSClimate <- climateData[Year>=1985 & Year<= 2011 & (Month > 9 | Month < 5),][
   ,.(NONGST = mean(Temp), NONGSP = sum(Prep), NONGSPET = sum(PET), NONGSCMI = sum(CMI), NONGSCO2 = mean(CO2)), 
   by = c("PlotID", "Year")]

allClimates <- setkey(AnnualClimate, PlotID, Year)[setkey(GSClimate, PlotID, Year), nomatch = 0]
allClimates <- setkey(allClimates, PlotID, Year)[setkey(NONGSClimate, PlotID, Year), nomatch = 0]

for(i in 1:nrow(inputData)){
  plotclimate <- allClimates[PlotID == inputData$PlotID[i],]
  iniyear <- inputData$IniYear[i]
  finyear <- inputData$FinYear[i]
  inputData$ATA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$AT)-mean(plotclimate$AT) 
  inputData$GSTA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$GST)-mean(plotclimate$GST)
  inputData$NONGSTA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$NONGST)-mean(plotclimate$NONGST)
  
  inputData$APA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$AP)-mean(plotclimate$AP) 
  inputData$GSPA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$GSP)-mean(plotclimate$GSP)
  inputData$NONGSPA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$NONGSP)-mean(plotclimate$NONGSP) 
  
  inputData$APETA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$APET)-mean(plotclimate$APET) 
  inputData$GSPETA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$GSPET)-mean(plotclimate$GSPET)
  inputData$NONGSPETA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$NONGSPET)-mean(plotclimate$NONGSPET) 
  
  inputData$ACMIA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$ACMI)-mean(plotclimate$ACMI) 
  inputData$GSCMIA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$GSCMI)-mean(plotclimate$GSCMI)
  inputData$NONGSCMIA[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$NONGSCMI)-mean(plotclimate$NONGSCMI)
  inputData$ACO2A[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$CO2)-mean(plotclimate$CO2) 
  inputData$GSCO2A[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$GSCO2)-mean(plotclimate$GSCO2)
  inputData$NONGSCO2A[i] <- mean(plotclimate[Year > iniyear & Year <= finyear,]$NONGSCO2)-mean(plotclimate$NONGSCO2)
  
}

write.csv(inputData, file.path(workPath, "data", "plotclimates.csv"), row.names=F)

workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "AllCensus_PositiveGrowth_RandomPlotADTree",
                                   "finalData.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table

thedata <- setkey(analysesData, PlotID, IniYear, FinYear)[setkey(inputData, PlotID, IniYear, FinYear),
                                                     nomatch = 0]

write.csv(thedata, file.path(workPath, "data", "newAllDataRescaledComp.csv"), row.names = FALSE)


