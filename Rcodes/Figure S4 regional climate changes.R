rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
temperature <- c("ATA", "GSTA", "NONGSTA")
allCMI <- c("ACMIA", "GSCMIA", "NONGSCMIA")
allP <- c("APA", "GSPA", "NONGSPA")
allPET <- c("APETA", "GSPETA", "NONGSPETA")
allCO2 <- c("ACO2A", "GSCO2A", "NONGSCO2A")
newlongcol <- c(temperature, allCMI, allP, allPET, allCO2)
climates <- read.csv(file.path(workPath, "data", "plotClimates.csv"), header = TRUE,
                     stringsAsFactors = FALSE) %>%
  data.table

climates <- climates[order(PlotID, IniYear),]
setnames(climates, "FinYear", "MidYear")
allclimatesvariables <- c("MidYear", temperature, allCMI, allP, allPET, allCO2)
allclimatesvariablesNew <- paste("Fin_", allclimatesvariables, sep = "")

climates[, (allclimatesvariablesNew) := shift(.SD, 1, NA, type = "lead"), 
         .SDcols = allclimatesvariables, by = PlotID]


climates <- climates[!is.na(Fin_MidYear), ]
for(i in 2:length(allclimatesvariables)){
  climates$tempV <- (climates[,allclimatesvariables[i], with = FALSE]+climates[,allclimatesvariablesNew[i], with = FALSE])/2
  set(climates, , c(allclimatesvariables[i], allclimatesvariablesNew[i]), NULL)
  setnames(climates, "tempV", allclimatesvariables[i])
}
setnames(climates, c("MidYear", "Fin_MidYear"), c("Year", "FinYear"))



selectionMethod <- "Year10Analyses"
analysesData <- fread(file.path(workPath, "data", selectionMethod, "finalData10.csv"))
analysesData <- analysesData[allCensusLiveTree == "yes",]
analysesDataPlotYear <- unique(analysesData[,.(PlotID, IniYear, FinYear)], by = c("PlotID", "IniYear"))

climates <- setkey(climates, PlotID, IniYear, FinYear)[setkey(analysesDataPlotYear, PlotID, IniYear, FinYear),
                                                       nomatch = 0]

climate_longform <- reshape(data = climates, varying = newlongcol, v.names = "Value",
                            times = newlongcol, timevar = "DependentVariable", 
                            direction = "long") %>% data.table
climate_longform[,Yearctd:=Year-mean(Year)]

for(indiclimategroup in c("temperature", "allCMI", "allP", "allPET", "allCO2")){
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
  if(indiclimategroup == "temperature"){ 
    allClimateData <- climateData
    allcoeff <- coeff
  } else {
    allClimateData <- rbind(allClimateData, climateData)
    allcoeff <- rbind(allcoeff, coeff)
  }
}

allClimateData$linetype <- factor(allClimateData$linetype, levels = c(1, 2))
allClimateData[Climate %in% c("ATA", "ACMIA", "APA", "APETA", "ACO2A"), Season:="Whole year"]
allClimateData[Climate %in% c("GSTA", "GSCMIA", "GSPA", "GSPETA", "GSCO2A"), Season:="Growing season"]
allClimateData[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSPA", "NONGSPETA","NONGSCO2A"), Season:="Non-growing season"]
allClimateData[, ':='(Season=factor(Season, levels = c("Whole year", "Growing season", "Non-growing season")),
                      ClimateName = factor(ClimateName, levels = c("temperature", "allCMI", "allP", "allPET", "allCO2")))]

relativePosition <- 0.1
slopestexts1 <- allClimateData[,.(maxValue = max(Value), minValue = min(Value)), by = ClimateName]
slopestexts1[, ':='(y = minValue+(maxValue-minValue)*relativePosition, Year = 1994, labels = "Slope:")]
allcoeff[DependentVariable %in% c("ATA", "ACMIA", "APA", "APETA", "ACO2A"), Season:="Whole year"]
allcoeff[DependentVariable %in% c("GSTA", "GSCMIA", "GSPA", "GSPETA", "GSCO2A"), Season:="Growing season"]
allcoeff[DependentVariable %in% c("NONGSTA", "NONGSCMIA", "NONGSPA", "NONGSPETA","NONGSCO2A"), Season:="Non-growing season"]
allcoeff[DependentVariable %in% temperature, ClimateName := "temperature"]
allcoeff[DependentVariable %in% allCMI, ClimateName := "allCMI"]
allcoeff[DependentVariable %in% allP, ClimateName := "allP"]
allcoeff[DependentVariable %in% allPET, ClimateName := "allPET"]
allcoeff[DependentVariable %in% allCO2, ClimateName := "allCO2"]
slopestexts1 <- left_join(allcoeff, slopestexts1[,.(ClimateName, y, Year)], by = "ClimateName") %>% data.table
slopestexts1[ClimateName == "temperature", labels := paste("slope:~", Value, "~degree~C~year^{-1}", sep = "")]
slopestexts1[ClimateName == "allCMI", labels := paste("slope:~", Value, "~mm~year^{-1}", sep = "")]
slopestexts1[ClimateName == "allP", labels := paste("slope:~", Value, "~mm~year^{-1}", sep = "")]
slopestexts1[ClimateName == "allPET", labels := paste("slope:~", Value, "~mm~year^{-1}", sep = "")]
slopestexts1[ClimateName == "allCO2", labels := paste("slope:~", Value, "~ppm~year^{-1}", sep = "")]
slopestexts1[, ':='(ClimateName = factor(ClimateName, levels = c("temperature", "allCMI", "allP", "allPET", "allCO2")),
                    Season=factor(Season, levels = c("Whole year", "Growing season", "Non-growing season")))]

anewlaberller <- list("temperature" = expression(atop("Temperature anomaly", paste("(", degree, "C)"))),
                      "allCMI" = expression(atop("CMI anomaly", "(mm)")),
                      "allP" = expression(atop("Precipitation anomaly", "(mm)")),
                      "allPET" = expression(atop("PET anomaly", "(mm)")),
                      'allCO2'= expression(atop(paste("C", O[2], " anomaly"), "(ppm)")))
anewlaberller2 <- list("Whole year" = "Annual",
                      "Growing season" = "Growing season", 
                      "Non-growing season" = "Non-growing season")
figurea_labeller <- function(variable,value){
  if(variable == "ClimateName"){
    return(anewlaberller[value])
  } else {
    return(anewlaberller2[value])
  }
}

FigureA <- ggplot(data = allClimateData[ClimateName != "allPET", ],
                  aes(x = Year, y = Value))+
  facet_grid(ClimateName~Season, switch = "y", scales = "free_y", labeller = figurea_labeller)+
  geom_point(col = "gray")+
  geom_line(aes(x = Year, y = predValue, linetype = linetype), size = 1)+
  scale_linetype(guide = "none")+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  geom_text(data = slopestexts1[ClimateName != "allPET", ], 
            aes(x = Year, y = y, label = labels), size = 5, parse = TRUE, hjust = 0)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 15, vjust = 2),
        strip.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 15, face = "italic"),
        strip.background = element_blank())
ggsave(file.path(workPath, "TablesFigures", "Figure S4. regional climate changes.png"), FigureA,
       width = 12, height = 10)

