rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(gridExtra)


workPath <- "~/GitHub/Climate_Growth"

selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"


analysesData <- fread(file.path(workPath, "data", selectionMethod, "finalData.csv"))
analysesData <- analysesData[allCensusLiveTree == "yes" & positiveGrowthTree == "yes",]
studySpecies <- c("All species", "Jack pine", "Trembling aspen",
                  "Black spruce", "Minor species")

for(indispecies in studySpecies){
  speciesdata <- analysesData[Species == indispecies,.(PlotID, uniTreeID, ABGR = BiomassGR,
                                                       logABGR = log(BiomassGR),
                                                         DBH = IniDBH, logDBH = log(IniDBH),
                                                        SA = IniFA, logSA = log(IniFA),
                                                        H, logH = log(H), 
                                                       Year, 
                                                         ATA, GSTA, NONGSTA,
                                                         APA, GSPA, NONGSPA,
                                                         ACMIA, GSCMIA, NONGSCMIA,
                                                         ACO2A, GSCO2A, NONGSCO2A)]
  tempcolnames <- names(speciesdata)[3:23]
  speciesdata <- reshape(data = speciesdata, varying = tempcolnames,
                         v.names = "Values",
                         timevar = "tempV",
                         direction = "long")
  speciesdata <- speciesdata[,.(mean = round(mean(Values), 2), sd = round(sd(Values), 2),
                                min = round(min(Values), 2),
                                max = round(max(Values), 2)), by = tempV]

  speciesdata[,Variable:=tempcolnames]
  logColNames <- tempcolnames[grep("log", tempcolnames)]
  speciesdata[!(Variable %in% logColNames), summary := paste(mean, " ± ", sd, "(",
                                                min, " to ", max, ")", sep = "")]
  speciesdata[Variable %in% logColNames, summary:=paste(round(exp(mean), 2), " ± ",
                                                        round(exp(sd), 2))]

  speciesdata <- speciesdata[,.(tempV, Variable, summary)]
  names(speciesdata)[3] <- indispecies
  if(indispecies == "All species"){
    Table1Output <- speciesdata
    rm(speciesdata)
  } else {
    Table1Output <- setkey(Table1Output, Variable, tempV)[setkey(speciesdata, Variable, tempV), nomatch = 0]
    rm(speciesdata)
  }
}
Table1Output <- Table1Output[order(tempV),]
Table1Output[,tempV:=NULL]
write.csv(Table1Output, file.path(workPath, "TablesFigures","table S1.csv"), row.names = FALSE)
