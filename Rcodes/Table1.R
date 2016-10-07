rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(gridExtra)

workPath <- "~/GitHub/Climate_Growth"

analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>%
  data.table

speciesdata <- data.table::copy(analysesData)
speciesdata <- speciesdata[,.(PlotID, uniTreeID, ABGR = BiomassGR, IniDBH, Year, Hegyi, 
                              RBI = RBI,
                              ATA, GSTA, NONGSTA,
                              APA, GSPA, NONGSPA,
                              ACMIA, GSCMIA, NONGSCMIA,
                              ACO2A, GSCO2A, NONGSCO2A)]
plotsummary <- data.table(Variable = c(-2, -1, 0), 
                          summary = c(length(unique(speciesdata$PlotID)),
                                      length(unique(speciesdata$uniTreeID)),
                                      nrow(speciesdata)))
tempcolnames <- names(speciesdata)[3:19]
speciesdata <- reshape(data = speciesdata, varying = tempcolnames,
                       v.names = "Values",
                       timevar = "Variable",
                       direction = "long")
speciesdata <- speciesdata[,.(mean = round(mean(Values), 2), sd = round(sd(Values), 2),
                              min = round(min(Values), 2),
                              max = round(max(Values), 2)), by = Variable]
speciesdata <- speciesdata[,.(Variable, 
                              summary = paste(mean, " ± ", sd, "(",
                                              min, " to ", max, ")", sep = ""))]
Table1Output <- rbind(speciesdata, plotsummary)
names(Table1Output)[2] <- "allData_summary"
rm(speciesdata)
studySpecies <- c("JP", "TA", "BS", "Other")

for(indispecies in studySpecies){
  if(indispecies == "Other"){
    speciesdata <- analysesData[!(Species %in% c("JP", "TA", "BS")),.(PlotID, uniTreeID, ABGR = BiomassGR,
                                                                      IniDBH, Year, Hegyi, 
                                                                      RBI = RBI,
                                                                      ATA, GSTA, NONGSTA,
                                                                      APA, GSPA, NONGSPA,
                                                                      ACMIA, GSCMIA, NONGSCMIA,
                                                                      ACO2A, GSCO2A, NONGSCO2A)]
  } else {
    speciesdata <- analysesData[Species == indispecies,.(PlotID, uniTreeID, ABGR = BiomassGR,
                                                         IniDBH, Year, Hegyi, 
                                                         RBI = RBI,
                                                         ATA, GSTA, NONGSTA,
                                                         APA, GSPA, NONGSPA,
                                                         ACMIA, GSCMIA, NONGSCMIA,
                                                         ACO2A, GSCO2A, NONGSCO2A)]
  }
  
  plotsummary <- data.table(Variable = c(-2, -1, 0), 
                            summary = c(length(unique(speciesdata$PlotID)),
                                        length(unique(speciesdata$uniTreeID)),
                                        nrow(speciesdata)))
  tempcolnames <- names(speciesdata)[3:19]
  speciesdata <- reshape(data = speciesdata, varying = tempcolnames,
                         v.names = "Values",
                         timevar = "Variable",
                         direction = "long")
  speciesdata <- speciesdata[,.(mean = round(mean(Values), 2), sd = round(sd(Values), 2),
                                min = round(min(Values), 2),
                                max = round(max(Values), 2)), by = Variable]
  speciesdata <- speciesdata[,.(Variable, 
                                summary = paste(mean, " ± ", sd, "(",
                                                min, " to ", max, ")", sep = ""))]
  speciesdata <- rbind(speciesdata, plotsummary)
  
  names(speciesdata)[2] <- paste(indispecies, "_", names(speciesdata)[2], sep = "")
  
  Table1Output <- setkey(Table1Output, Variable)[setkey(speciesdata, Variable), nomatch = 0]
  
}
Table1Output <- Table1Output[order(Variable),]
Table1Output$Variable <- c("NofPlot", "NofTree", "NofObs", tempcolnames)
write.csv(Table1Output, file.path(workPath, "TablesFigures","table1.csv"), row.names = FALSE)
