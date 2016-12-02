rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(gridExtra)

workPath <- "~/GitHub/Climate_Growth"

analysesData <- read.csv(file.path(workPath, "data", "newAllDataRescaledComp.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>%
  data.table

studySpecies <- c("All species", "Jack pine", "Trembling aspen",
                  "Black spruce", "Other species")

for(indispecies in studySpecies){
  speciesdata <- analysesData[DataType == indispecies,.(PlotID, uniTreeID, ABGR = BiomassGR,
                                                        logABGR = log(BiomassGR), 
                                                         DBH = IniDBH, logDBH = log(IniDBH),
                                                        SA, logSA = log(SA),
                                                        Year, IntraH, logIntraH = log(IntraH+1),
                                                        InterH, logInterH = log(InterH+1), 
                                                         ATA, GSTA, NONGSTA,
                                                         APA, GSPA, NONGSPA,
                                                         ACMIA, GSCMIA, NONGSCMIA,
                                                         ACO2A, GSCO2A, NONGSCO2A)]
  tempcolnames <- names(speciesdata)[3:25]
  speciesdata <- reshape(data = speciesdata, varying = tempcolnames,
                         v.names = "Values",
                         timevar = "tempV",
                         direction = "long")
  speciesdata <- speciesdata[,.(mean = round(mean(Values), 2), sd = round(sd(Values), 2),
                                min = round(min(Values), 2),
                                max = round(max(Values), 2)), by = tempV]

  speciesdata[,Variable:=tempcolnames]
  speciesdata[, summary := paste(mean, " ± ", sd, "(",
                                                min, " to ", max, ")", sep = "")]
  speciesdata[Variable %in% c("logABGR", "logDBH", "logSA"), summary:=paste(round(exp(mean), 2), 
                                                                            "(", paste(round(exp(mean-sd), 2)),
                                                                            " to ", paste(round(exp(mean+sd), 2)),
                                                                            ")", sep = "")]
  speciesdata[Variable %in% c("logIntraH", "logInterH"), summary:=paste(round(exp(mean)-1, 2), 
                                                                            "(", paste(round(exp(mean-sd)-1, 2)),
                                                                            " to ", paste(round(exp(mean+sd)-1, 2)),
                                                                            ")", sep = "")]
  
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
