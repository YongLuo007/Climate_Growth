# to produce relative importance figure
rm(list = ls())
library(ggplot2, lib.loc = "~/GitHub/Climate_Growth/RequiredRPackages")
library(SpaDES)
library(data.table)



workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "Year10Analyses"
load(file.path(workPath, "data",selectionMethod,
               "FullYearModels.RData"))
allCoeffs <- data.frame(Variable = c("(Intercept)", "DBH", "H",
                                     "SA", "Year", 
                                     "DBH:H", "DBH:SA",
                                     "H:SA", "DBH:Year", 
                                     "Year:H", "Year:SA"),
                        stringsAsFactors = FALSE)

##### for overall temporal trends and its dependency on DBH and RBI
for(i in 1:length(fixedCoeffAll)){
  indicoeff <- fixedCoeffAll[[i]]
  indicoeff[, Variable:= unlist(lapply(lapply(lapply((rn), function(x) unlist(strsplit(x, ":"))), function(y) sort(y)),
                                       function(z) paste(z, collapse = ":")))]
  indicoeff[, Variable:= gsub("ctd", "", rn)]
  indicoeff[, Variable:= gsub("log", "", Variable)]
  indicoeff <- indicoeff[,.(Variable, EstimateSE = paste(round(Value, 4), "(",
                                                         round(`Std.Error`, 4), ")", sep = ""), 
                            P = round(`p-value`, 2))]
  names(indicoeff)[2:3] <- paste(names(fixedCoeffAll)[i], names(indicoeff)[2:3])
  if(i == 1){
    allcoeffs <- indicoeff
  } else {
    allcoeffs <- left_join(allcoeffs, indicoeff, by = "Variable")
  }
}

write.csv(allcoeffs, file.path(workPath, "TablesFigures", "Table1.csv"),
          row.names = FALSE)
