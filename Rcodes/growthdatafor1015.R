
rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table);library(MuMIn)
library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
analysesDataOrig <- fread(file.path(workPath, "data", selectionMethod, "finalData.csv"))
analysesDataOrig <- analysesDataOrig[Species != "All species",]
analysesDataOrig <- analysesDataOrig[order(uniTreeID, Year),]
analysesDataOrig[ ,':='(MidYear = FinYear, 
                        MidDBH = FinDBH, 
                        MidFA = FinFA,
                        MidBiomass = FinBiomass,
                        MidBA = FinBA)]
analysesDataOrig <- analysesDataOrig[,c("FinYear2", "FinDBH2", "FinFA2",  
                                        "FinBiomass2", "tempuniTreeID",
                                        "ATA2", "GSTA2", "NONGSTA2",
                                        "ACMIA2", "GSCMIA2",
                                        "ACO2A2", "FinBA2", "MidH") 
                                     := data.table::shift(x = analysesDataOrig[,.(FinYear, FinDBH, FinFA,
                                                                                  FinBiomass, uniTreeID,
                                                                                  ATA, GSTA, NONGSTA,
                                                                                  ACMIA, GSCMIA,
                                                                                  ACO2A, FinBA, H)],
                                                          n = 1, fill = NA, type = "lead", give.names = FALSE)]
analysesDataOrig <- analysesDataOrig[FinYear != FinYear2,]
analysesDataOrig <- analysesDataOrig[tempuniTreeID == uniTreeID,][,tempuniTreeID := NULL]
analysesDataOrig[, plotcensuslength:=length(unique(FinYear2)), by = PlotID]
analysesDataOrig <- analysesDataOrig[plotcensuslength>1]
analysesDataOrig[, ':='(BiomassGR2=(FinBiomass2-IniBiomass)/10,
                        BAGR2 = (FinBA2-IniBA)/10)]
analysesDataOrig <- rbind(data.table::copy(analysesDataOrig),
                          data.table::copy(analysesDataOrig)[,Species:="All species"])
analysesDataOrig[,':='(FinYear = FinYear2,
                       FinDBH = FinDBH2,
                       FinFA = FinFA2,
                       BiomassGR = BiomassGR2,
                       FinBA = FinBA2,
                       BAGR = BAGR2,
                       Year = (IniYear+FinYear2)/2,
                       FinBiomass = FinBiomass2,
                       ATA = (ATA+ATA2)/2,
                       GSTA = (GSTA+GSTA2)/2,
                       NONGSTA = (NONGSTA+NONGSTA2)/2,
                       ACMIA = (ACMIA+ACMIA2)/2,
                       GSCMIA = (GSCMIA+GSCMIA2)/2,
                       ACO2A = (ACO2A+ACO2A2)/2)]
set(analysesDataOrig, , c("FinYear2", "FinDBH2", "FinFA2",  
                          "FinBiomass2", 
                          "ATA2", "GSTA2", "NONGSTA2",
                          "ACMIA2", "GSCMIA2",
                          "ACO2A2", "APA", "GSPA",
                          "NONGSPA", "APETA", "GSPETA", 
                          "NONGSPETA", "NONGSCMIA",
                          "GSCO2A", "NONGSCO2A", "positiveGrowthTree",
                          "FinBA2", "plotcensuslength",
                          "BiomassGR2", "Year", "BAGR2"), NULL)
analysesDataOrig <- analysesDataOrig[,.(PlotID, uniTreeID, allCensusLiveTree, Species,
                                        species, IniYear, IniFA, IniDBH, IniH = H, IniBiomass, IniBA,
                                        MidYear, MidFA, MidDBH, MidH, MidBiomass, MidBA,
                                        FinYear, FinFA, FinDBH, FinBiomass, FinBA,
                                        BiomassGR, BAGR,
                                        ATA, GSTA, NONGSTA, ACMIA, GSCMIA, ACO2A)]
analysesDataOrig[IniDBH == FinDBH, ':='(BiomassGR = 0,
                                        BAGR = 0)]
analysesDataOrig <- analysesDataOrig[order(Species, PlotID, uniTreeID, IniYear),]


write.csv(analysesDataOrig, 
          file.path(workPath, "data", 
                    "Year10Analyses", "finalData10.csv"),
          row.names = F)



rm(list = ls())
library(dplyr); library(SpaDES); library(nlme); library(data.table);library(MuMIn)
library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
analysesDataOrig <- fread(file.path(workPath, "data",
                                    selectionMethod, "finalData10.csv"))
analysesDataOrig <- analysesDataOrig[Species != "All species",]
analysesDataOrig <- analysesDataOrig[order(uniTreeID, Year),]
analysesDataOrig <- analysesDataOrig[,c("FinYear2", "FinDBH2", "FinFA2",  
                                        "FinBiomass2", "tempuniTreeID",
                                        "ATA2", "GSTA2", "NONGSTA2",
                                        "ACMIA2", "GSCMIA2",
                                        "ACO2A2", "FinBA2") 
                                     := data.table::shift(x = analysesDataOrig[,.(FinYear, FinDBH, FinFA,
                                                                                  FinBiomass, uniTreeID,
                                                                                  ATA, GSTA, NONGSTA,
                                                                                  ACMIA, GSCMIA,
                                                                                  ACO2A, FinBA)],
                                                          n = 1, fill = NA, type = "lead", give.names = FALSE)]
analysesDataOrig <- analysesDataOrig[FinYear != FinYear2,]
analysesDataOrig <- analysesDataOrig[tempuniTreeID == uniTreeID,][,tempuniTreeID := NULL]
analysesDataOrig[, plotcensuslength:=length(unique(FinYear2)), by = PlotID]
analysesDataOrig <- analysesDataOrig[plotcensuslength>1]
analysesDataOrig[, BiomassGR2:=(FinBiomass2-IniBiomass)/10]
analysesDataOrig <- rbind(data.table::copy(analysesDataOrig),
                          data.table::copy(analysesDataOrig)[,Species:="All species"])
analysesDataOrig[,':='(FinYear = FinYear2,
                       FinDBH = FinDBH2,
                       FinFA = FinFA2,
                       BiomassGR = BiomassGR2,
                       FinBA = FinBA2,
                       Year = (IniYear+FinYear2)/2,
                       FinBiomass = FinBiomass2,
                       ATA = (ATA+ATA2)/3,
                       GSTA = (GSTA+GSTA2)/3,
                       NONGSTA = (NONGSTA+NONGSTA2)/2,
                       ACMIA = (ACMIA+ACMIA2)/2,
                       GSCMIA = (GSCMIA+GSCMIA2)/2,
                       ACO2A = (ACO2A+ACO2A2)/2)]
set(analysesDataOrig, , c("FinYear2", "FinDBH2", "FinFA2",  
                          "FinBiomass2", 
                          "ATA2", "GSTA2", "NONGSTA2",
                          "ACMIA2", "GSCMIA2",
                          "ACO2A2", "APA", "GSPA",
                          "NONGSPA", "APETA", "GSPETA", 
                          "NONGSPETA", "NONGSCMIA",
                          "GSCO2A", "NONGSCO2A", "positiveGrowthTree",
                          "FinBA2", "plotcensuslength",
                          "BiomassGR2"), NULL)
write.csv(analysesDataOrig, file.path(workPath, "data", 
                                      selectionMethod, "finalData15.csv"),
          row.names = F)


