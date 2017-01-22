rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- "J:/MBgrowth"
}
petMethods <- c("pen", "pt", "rfc", "rfu", "rg")
for(petmethod in petMethods){
  inputData <- read.csv(file.path(workPath, "data", "StudyAreaClimates_BiomSIM", "Loc", "selectedplots.csv"), header = TRUE,
                        stringsAsFactors = FALSE) %>%
    data.table
  petData <- read.csv(file.path(workPath, "data", "newClimatesRobbie", paste("NACID_ForLuo_ManitobaPlots_etp_", 
                                                                             petmethod, ".csv", sep = "")),
                      header = FALSE, stringsAsFactors = FALSE) %>% data.table
  names(petData) <- inputData$Name
  prepData <- read.csv(file.path(workPath, "data", "newClimatesRobbie", "NACID_ForLuo_ManitobaPlots_prcp.csv"),
                       header = FALSE, stringsAsFactors = FALSE) %>% data.table
  names(prepData) <- inputData$Name
  timeData <- read.csv(file.path(workPath, "data", "newClimatesRobbie", "NACID_ForLuo_ManitobaPlots_time.csv"),
                       header = FALSE, stringsAsFactors = FALSE) %>% data.table
  names(timeData) <- c("Year", "Month")
  timeData$NofDay <- as.numeric(diff(seq(as.Date("1980-01-01"), as.Date("2015-01-01"), by = "month")))
  petData <- cbind(timeData, petData)
  
  petDataLong <- melt(petData, id.vars = c("Year", "Month", "NofDay"),
                      measure.vars = inputData$Name, variable.name = "PlotID",
                      value.name = "PET")
  petDataLong[, PET:=NofDay*PET]
  petDataLong[, NofDay:=NULL]
  prepData <- cbind(timeData, prepData)
  prepDataLong <- melt(prepData, id.vars = c("Year", "Month", "NofDay"),
                       measure.vars = inputData$Name, variable.name = "PlotID",
                       value.name = "Prep")
  prepDataLong[, Prep:=NofDay*Prep]
  prepDataLong[, NofDay:=NULL]
  climateall <- setkey(petDataLong, Year, Month, PlotID)[setkey(prepDataLong, Year, Month, PlotID),
                                                         nomatch = 0]
  CO2data <- read.csv(file.path(workPath, "data", "StudyAreaClimates_BiomSIM", "CO2.csv"),
                      header = TRUE, stringsAsFactors = FALSE) %>% data.table
  
  climateData <- setkey(climateall, Year, Month)[setkey(CO2data[,.(Year, Month, CO2)], Year, Month),
                                                 nomatch = 0]
  climateData[,Year:=as.numeric(Year)]
  climateData[,CMI:=Prep - PET]
  climateData[Month>=10, Year:=as.integer(Year)+1]
  AnnualClimate <- climateData[Year>=1984 & Year<= 2011,][,.(AP = sum(Prep),
                                                             APET = sum(PET),
                                                             ACMI = sum(CMI), CO2 = mean(CO2)), 
                                                          by = c("PlotID", "Year")]
  GSClimate <- climateData[Year>=1984 & Year<= 2011 & Month<=9 & Month >= 5,][
    ,.(GSP = sum(Prep), GSPET = sum(PET), GSCMI = sum(CMI), GSCO2 = mean(CO2)), by = c("PlotID", "Year")]
  
  NONGSClimate <- climateData[Year>=1984 & Year<= 2011 & (Month > 9 | Month < 5),][
    ,.(NONGSP = sum(Prep), NONGSPET = sum(PET), NONGSCMI = sum(CMI), NONGSCO2 = mean(CO2)), 
    by = c("PlotID", "Year")]
  
  allClimates <- setkey(AnnualClimate, PlotID, Year)[setkey(GSClimate, PlotID, Year), nomatch = 0]
  allClimates <- setkey(allClimates, PlotID, Year)[setkey(NONGSClimate, PlotID, Year), nomatch = 0]
  
  inputData <- read.csv(file.path(workPath, "data", "plotsummary.csv"), header = TRUE, 
                        stringsAsFactors = FALSE) %>% data.table
  
  for(i in 1:nrow(inputData)){
    plotclimate <- allClimates[PlotID == inputData$PlotID[i],]
    iniyear <- inputData$IniYear[i]
    finyear <- inputData$FinYear[i]
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


  allCMI <- c("ACMIA", "GSCMIA", "NONGSCMIA")
  allP <- c("APA", "GSPA", "NONGSPA")
  allPET <- c("APETA", "GSPETA", "NONGSPETA")
  allCO2 <- c("ACO2A", "GSCO2A", "NONGSCO2A")
  newlongcol <- c(allCMI, allP, allPET, allCO2)
  biosim <- TRUE
  if(biosim){
    climates <- read.csv(file.path(workPath, "data", "plotclimates.csv"),
                         header = TRUE, stringsAsFactors = FALSE) %>% data.table
  } else {
    climates <- data.table::copy(inputData)
  }
  
  climates[,Year:=(FinYear+IniYear)/2]
  climate_longform <- reshape(data = climates, varying = newlongcol, v.names = "Value",
                              times = newlongcol, timevar = "DependentVariable", 
                              direction = "long") %>% data.table
  climate_longform[,Yearctd:=Year-mean(Year)]
  
  for(indiclimategroup in c("allCMI", "allP", "allPET", "allCO2")){
    climateData <- climate_longform[DependentVariable %in% get(indiclimategroup),]
    climateModel <- lme(Value ~ DependentVariable/Yearctd, random =~(DependentVariable-1)|PlotID,
                        data = climateData,
                        control = lmeControl(opt="optim", maxIter=50000, msMaxIter = 50000))
    climateData$predValue <- predict(climateModel, newdata = climateData, level = 0) 
    coeff <- data.frame(summary(climateModel)$tTable)
    coeff <- coeff[row.names(coeff) %in% paste("DependentVariable", get(indiclimategroup), ":Yearctd", sep = ""),c(1,5)] %>%
      data.table
    coeff[,':='(DependentVariable = get(indiclimategroup),
                Value = round(coeff$Value, 2), linetype=1)]
    coeff[p.value>=0.05, linetype:=2]
    climateData <- setkey(climateData, DependentVariable)[setkey(coeff[,.(DependentVariable, linetype)], DependentVariable),
                                                          nomatch = 0]
    climateData <- climateData[,.(ClimateName = indiclimategroup, Climate = DependentVariable,
                                  Year, Value, predValue, linetype)]
    if(indiclimategroup == "allCMI"){ 
      allClimateData <- climateData
      allcoeff <- coeff[,petMethod:=petmethod]
    } else {
      allClimateData <- rbind(allClimateData, climateData)
      allcoeff <- rbind(allcoeff, coeff[,petMethod:=petmethod])
    }
  }
  if(petmethod == "pen"){
    allmethodcoeff <- allcoeff
  } else {
    allmethodcoeff <- rbind(allmethodcoeff, allcoeff)
  }
}
write.csv(allmethodcoeff, file.path(workPath, "data", "pettrends.csv"),row.names = F)



