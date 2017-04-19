# competition relative to assymetric and crowdedness
rm(list = ls())
library(SpaDES); library(data.table)

workPath <- "~/Github/Climate_Growth"
library(ggplot2, lib.loc = file.path(workPath, "RequiredRPackages"))
competitionDataOrig <- fread((file.path(workPath, "data", "Year10Analyses", "forcompetitionIndex.csv")))



# "AS" 
competitionDataOrig[Species == "AS",Species_Std:="black spruce"]
# "BA"
competitionDataOrig[Species == "BA",Species_Std:="balsam poplar"]
# "BF" 
competitionDataOrig[Species == "BF",Species_Std:="balsam fir"]
# "BO"
competitionDataOrig[Species == "BO",Species_Std:="white oak"] # bur oak to
# "BS"
competitionDataOrig[Species == "BS",Species_Std:="black spruce"]
# "EC"
competitionDataOrig[Species == "EC",Species_Std:="western redcedar"] # cedar
# "JP"
competitionDataOrig[Species == "JP",Species_Std:="jack pine"]
# "MM"
competitionDataOrig[Species == "MM",Species_Std:="silver maple"]
# "RP"
competitionDataOrig[Species == "RP",Species_Std:="red pine"]
# "TA"
competitionDataOrig[Species == "TA",Species_Std:="trembling aspen"]
# "TL"
competitionDataOrig[Species == "TL",Species_Std:="tamarack larch"]
# "WB"
competitionDataOrig[Species == "WB",Species_Std:="white birch"]
# "WE"
competitionDataOrig[Species == "WE",Species_Std:="white elm"]
# "WS"
competitionDataOrig[Species == "WS",Species_Std:="white spruce"]

source(file.path(workPath, "Rcodes",  "Rfunctions", "biomassCalculation.R"))

competitionDataOrig$Biomass <- biomassCalculation(species = competitionDataOrig$Species_Std,
                                         DBH = competitionDataOrig$DBH)

competitionDataOrig[, ':='(minBiomass = min(Biomass),
                           maxBiomass = max(Biomass)),
                    by = c("PlotID", "Year")]
competitionData <- competitionDataOrig[,.(relaBiomass = (maxBiomass-Biomass)/(maxBiomass-minBiomass)), 
                                       by = c("PlotID", "Year", "uniTreeID")]
setnames(competitionData, "Year", "MidYear")

analysesData <- fread(file.path(workPath, "data", "Year10Analyses", "finalData10.csv"))
analysesData <- analysesData[allCensusLiveTree == "yes",]
studySpeciesGroup <- c("All species", "Jack pine", "Trembling aspen",
                       "Black spruce", "Minor species")
analysesDataH <- analysesData[,.(Species_Group, PlotID, uniTreeID, MidBiomass, MidYear, H)]

for(indispecies in studySpeciesGroup){
  indispdata <- analysesDataH[Species_Group == indispecies,]
  indispdata <- setkey(indispdata, PlotID, uniTreeID, MidYear)[setkey(competitionData, PlotID,
                                                                      uniTreeID, MidYear),
                                                               nomatch = 0]
  if(indispecies == "All species" | indispecies ==  "Jack pine"){
    j <- 1.2
  } else if (indispecies == "Trembling aspen"){ 
    j <- 0.9
  } else if (indispecies == "Black spruce"){
    j <- 1.3 
  } else {
    j <- 0.5
  }
  crowdH <- fread(file.path(workPath, "data", "AllCompetitionData",
                               paste("CompetitionData_DW",j, "_SW", 0, ".csv", sep = "")))
  crowdH <- crowdH[,.(uniTreeID, MidYear = IniYear, crowdH = H)]
  indispdata <- setkey(indispdata, uniTreeID, MidYear)[setkey(crowdH, uniTreeID, MidYear),
                                                       nomatch = 0]
  indispdata[, crowH := MidBiomass*crowdH]
  indispdata[, assymC := H/crowH]
  if(indispecies == "All species"){
    alldata <- indispdata
  } else {
    alldata <- rbind(alldata, indispdata)
  }
}
figuredata <- alldata[,.(Species_Group, uniTreeID, MidYear, H, crowdC = crowH, assymC)]

figuredata <- melt(figuredata, id.vars = c("Species_Group", "uniTreeID", "MidYear", "H"),
                   measure.vars = c("crowdC", "assymC"))
coroutput <- data.table::data.table(Species_Group = character(),
                                    CompetitionType = character(),
                                    Cor = numeric())
for(indispecies in studySpeciesGroup){
  figuredataindi <- figuredata[Species_Group == indispecies,]
  for(j in c("crowdC", "assymC")){
    figuredataindi_indic <- figuredataindi[variable == j,]
    coroutputindi <- data.table(Species_Group = indispecies, 
                                CompetitionType = j,
                                Cor = cor(figuredataindi_indic$H,
                                          figuredataindi_indic$value))
    coroutput <- rbind(coroutput, coroutputindi)
  }
}
coroutput[,':='(Species_Group = factor(Species_Group, levels = studySpeciesGroup),
                CompetitionType = factor(CompetitionType, levels = c("crowdC", "assymC"),
                                         labels = c("log(Crowdedness competition)",
                                                    "log(Asymmetric competition)")),
                text = paste("Correlation coefficient:", round(Cor, 2)),
                x = -Inf, y = -Inf)]

figuredata[,':='(Species_Group = factor(Species_Group, levels = studySpeciesGroup),
                 CompetitionType = factor(variable, levels = c("crowdC", "assymC"),
                                          labels = c("log(Crowdedness competition)",
                                                     "log(Asymmetric competition)")))]
texttable <- data.table(Species_Group = factor(studySpeciesGroup, levels = studySpeciesGroup),
                        CompetitionType = factor("log(Crowdedness competition)",
                                                 levels = c("log(Crowdedness competition)",
                                                            "log(Asymmetric competition)")),
                        value = -Inf, H = Inf, labels = letters[1:5])
segtable <- rbind(data.table(Species_Group = factor(studySpeciesGroup, levels = studySpeciesGroup),
                             CompetitionType = factor("log(Crowdedness competition)",
                                                      levels = c("log(Crowdedness competition)",
                                                                 "log(Asymmetric competition)")),
                             x = -Inf, xend = -Inf, y = -Inf, yend = Inf),
                  data.table(expand.grid(Species_Group = factor(studySpeciesGroup, 
                                                                levels = studySpeciesGroup),
                                         CompetitionType = factor(c("log(Crowdedness competition)",
                                                                    "log(Asymmetric competition)"),
                                                                  levels = c("log(Crowdedness competition)",
                                                                             "log(Asymmetric competition)"))))[, ':='(x = -Inf, xend = Inf, y = -Inf, yend = -Inf)])




thefigure <- ggplot(data = figuredata, aes(x = log(value), y = log(H)))+
  facet_grid(Species_Group~CompetitionType, scales = "free_x", switch = "x")+
  geom_point()+
  geom_text(data = texttable, aes(x = value, y = H, label = labels), 
            size = 7, hjust = -0.5, vjust = 1)+
  geom_text(data = coroutput, aes(x = x, y = y, label = text),
            size = 4, hjust = -0.35, vjust = -0.5, col = "red")+
  geom_segment(data = segtable, aes(x = x, xend = xend, y = y, yend = yend))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_blank())

ggsave(file.path(workPath, "TablesFigures", "CompetitionRelativeness.png"),
       thefigure, width = 7, height = 9)
  










