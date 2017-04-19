rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(gridExtra)


workPath <- "~/GitHub/Climate_Growth"

selectionMethod <- "Year10Analyses"


analysesData <- fread(file.path(workPath, "data", selectionMethod, "finalData10.csv"))
analysesData <- analysesData[allCensusLiveTree == "yes",]
studySpecies <- c("All species", "Jack pine", "Trembling aspen",
                  "Black spruce", "Minor species")

for(indispecies in studySpecies){
  speciesdata <- analysesData[Species_Group == indispecies,.(PlotID, uniTreeID, ABGR = BiomassGR,
                                                       DBH = MidDBH, logDBH = log(MidDBH),
                                                       SA = MidFA, logSA = log(MidFA),
                                                       H = H, logH = log(H), 

                                                       Year = MidYear, 
                                                       ATA, GSTA, NONGSTA,
                                                       ACMIA, GSCMIA, 
                                                       ACO2A)]

  tempcolnames <- names(speciesdata)[3:length(names(speciesdata))]
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
                                                min, " ~ ", max, ")", sep = "")]
  speciesdata[Variable %in% logColNames, summary:=paste(round(exp(mean), 2), " ± ",
                                                        round(exp(sd), 2))]

  speciesdata <- speciesdata[,.(tempV, Variable, summary)]
  names(speciesdata)[3] <- indispecies
  speciesdataaddon <- analysesData[Species_Group == indispecies,
                                   .(tempV = length(tempcolnames)+1,
                                     Variable = "Q95H",
                                     Q95H = paste(round(quantile(H, 0.025), 2),
                                                  round(quantile(H, 0.975), 2),
                                                  sep = " ~ "))]
  
  names(speciesdataaddon)[3] <- indispecies
  speciesdata <- rbind(speciesdata, speciesdataaddon)
  
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
