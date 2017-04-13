rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
thedata <- read.csv(file.path(workPath, "data", "Year10Analyses",
                              "MBdataSimplified.csv"), header = TRUE,
                    stringsAsFactors = FALSE) %>%
  data.table
source(file.path(workPath, "Rcodes",  "Rfunctions", "biomassCalculation.R"))
thedata$IniBiomass <- biomassCalculation(species = thedata$Species_Std,
                                         DBH = thedata$IniDBH)
thedata$MidBiomass <- biomassCalculation(species = thedata$Species_Std,
                                         DBH = thedata$MidDBH)
thedata$FinBiomass <- biomassCalculation(species = thedata$Species_Std,
                                         DBH = thedata$FinDBH)
thedata[,':='(IniBA = 3.1415926*((IniDBH/2)^2),
              FinBA = 3.1415926*((FinDBH/2)^2))]
thedata[,':='(BAGR = (FinBA-IniBA)/(FinYear-IniYear),
              BiomassGR = (FinBiomass - IniBiomass)/(FinYear-IniYear))]
thedata[, Species_Group := "Minor species"]
thedata[Species_Std == "jack pine", Species_Group := "Jack pine"]
thedata[Species_Std == "trembling aspen", Species_Group := "Trembling aspen"]
thedata[Species_Std == "black spruce", Species_Group := "Black spruce"]


write.csv(thedata,
          file.path(workPath, "data", "Year10Analyses", "MBdatafinal.csv"),
          row.names = FALSE)


