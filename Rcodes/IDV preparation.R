rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
thedata <- read.csv(file.path(workPath, "data", "MBdataSimplified.csv"), header = TRUE,
                    stringsAsFactors = FALSE) %>%
  data.table
source(file.path(workPath, "Rcodes",  "Rfunctions", "biomassCalculation.R"))
thedata$IniBiomass <- biomassCalculation(species = thedata$species,
                                         DBH = thedata$IniDBH)
thedata$FinBiomass <- biomassCalculation(species = thedata$species,
                                         DBH = thedata$FinDBH)
thedata[,':='(IniBA = 3.1415926*((IniDBH/2)^2), 
              FinBA = 3.1415926*((FinDBH/2)^2))]
thedata[,':='(BAGR = (FinBA-IniBA)/(FinYear-IniYear),
              BiomassGR = (FinBiomass - IniBiomass)/(FinYear-IniYear),
              Year = (FinYear+IniYear)/2)]
thedata[Species == "JP", DataType := "Jack pine"]
thedata[Species == "TA", DataType := "Trembling aspen"]
thedata[Species == "BS", DataType := "Black spruce"]
thedata[!(Species %in% c("JP", "TA", "BS")), DataType := "Other species"]
thedata[,Species:=DataType]
thedata[,DataType:=NULL]
write.csv(thedata,
          file.path(workPath, "data", "MBdatafinal.csv"),
          row.names = FALSE)


