rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(parallel)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Github/Climate_Growth"
} else {
  workPath <- file.path("", "home", "yonluo","Climate_Growth")
}
firstRun <- FALSE
if(firstRun){
  load(file.path(workPath, "data", "MBPSP.RData"))
  names(allPSP) <- c("PlotID", "FMU", "TWP", "RGE", "PlotNumber",
                     "StandType", "StructureType", "VegeType",
                     "SStructure", "Moist", "YearEstabish", "Easting",
                     "Northing", "TreeNumber", "Species", "Distance",
                     "Angle", "DBH", "Height", "Status", "Class",
                     "TreeAge", "Year", "PlotSize")
  inputData <- data.table::copy(allPSP)
  analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                           stringsAsFactors = FALSE) %>% data.table
  inputData <- inputData[PlotID %in% unique(analysesData$PlotID),]
  
  inputData[Species == "AS",species:="black spruce"]
  # "BA"
  inputData[Species == "BA",species:="balsam poplar"]
  # "BF" 
  inputData[Species == "BF",species:="balsam fir"]
  # "BO"
  inputData[Species == "BO",species:="white oak"] # bur oak to
  # "BS"
  inputData[Species == "BS",species:="black spruce"]
  # "EC"
  inputData[Species == "EC",species:="western redcedar"] # cedar
  # "JP"
  inputData[Species == "JP",species:="jack pine"]
  # "MM"
  inputData[Species == "MM",species:="silver maple"]
  # "RP"
  inputData[Species == "RP",species:="red pine"]
  # "TA"
  inputData[Species == "TA",species:="trembling aspen"]
  # "TL"
  inputData[Species == "TL",species:="tamarack larch"]
  # "WB"
  inputData[Species == "WB",species:="white birch"]
  # "WE"
  inputData[Species == "WE",species:="white elm"]
  # "WS"
  inputData[Species == "WS",species:="white spruce"]
  # for unknownspecies
  inputData[is.na(species), species:="unknown"]
  
  inputData <- inputData[Status == 1 | Status == 0,]
  source(file.path(workPath, "Rcodes",  "Rfunctions", "biomassCalculation.R"))
  inputData$Biomass <- biomassCalculation(species = inputData$species,
                                          DBH = inputData$DBH)
  source(file.path(workPath, "Rcodes", "Rfunctions", "HeghyiCICalculationModified.R"))
  inputData[PlotID == "46-132", Distance:=Distance/100]
  CIcompetitionData <- inputData[!is.na(Distance),]
  CIcompetitionData <- CIcompetitionData[,.(PlotNumber = PlotID, TreeNumber = TreeNumber,
                                            Year, Distance, Angle, DBH, Species, Biomass)]
  set(analysesData, , c("H_DBH", "IntraH_DBH", "InterH_DBH", 
                        "H_Biomass", "IntraH_Biomass", "InterH_Biomass", "IntraH1_3", 
                        "InterH0_4", "H"), NULL)
  
  i <-  c(0.5, 0.6, 0.7, 1, 1.1)
  j  <-  c(0.4, 0.7, 0.8, 1, 4, 5.3)
  processCIdata <- data.table::copy(CIcompetitionData)
  CIdata <- HeghyiCICalculation(data = processCIdata,
                                maxRadius = 12.62,
                                sizeIndex = "Biomass",
                                distanceWeight = i, 
                                sizeWeight = j,
                                assymetricScale = "Rescale",
                                testing = TRUE)
  verifydata <- CIdata$verifyOutput
  verifydata[,':='(uniTreeID = paste(PlotNumber, "_", TreeNumber, sep = ""))]
  verifydata <- unique(verifydata, by = c("uniTreeID", "Year"))
  analysesData <- read.csv(file.path(workPath, "data", "newAllDataRescaledComp.csv"), header = TRUE,
                           stringsAsFactors = FALSE) %>% data.table
  
  for(indispecies in c("Jack pine", "Trembling aspen", "Black spruce")){
    for(competitionName in c("Intraspecific competition", "Interspecific competition")){
      speciesdata <- analysesData[DataType == indispecies,.(uniTreeID, Year = IniYear)]
      if(indispecies == "Jack pine" & competitionName == "Intraspecific competition"){
        a <- 0.4
        b <- 1.1
      } else if (indispecies == "Jack pine" & competitionName == "Interspecific competition"){
        a <- 5.3
        b <- 1
      } else if (indispecies == "Trembling aspen" & competitionName == "Intraspecific competition"){
        a <- 0.8
        b <- 0.5
      } else if (indispecies == "Trembling aspen" & competitionName == "Interspecific competition"){
        a <- 4
        b <- 0.7
      } else if (indispecies == "Black spruce" & competitionName == "Intraspecific competition"){
        a <- 1
        b <- 1
      } else if (indispecies == "Black spruce" & competitionName == "Interspecific competition"){
        a <- 0.7
        b <- 0.6
      }
      if(competitionName == "Intraspecific competition"){
        HData <- verifydata[, c("uniTreeID", "Year", "areaRatio", paste(c("correctedIntraH_", "obsIntraH_"),
                                                                        "DW", b, "_SW",a, 
                                                                        sep = "")), with = FALSE]
      } 
      if(competitionName == "Interspecific competition") {
        HData <- verifydata[, c("uniTreeID", "Year", "areaRatio", paste(c("correctedInterH_", "obsInterH_"),
                                                                        "DW", b, "_SW",a, 
                                                                        sep = "")), with = FALSE]
      }
      names(HData)[4:5] <- c("corrected", "observed")
      speciesdata <- setkey(speciesdata, uniTreeID, Year)[setkey(HData, uniTreeID, Year),
                                                          nomatch = 0]
      rm(HData)
      speciesdata <- speciesdata[,.(Species = indispecies, uniTreeID, Year, areaRatio, CompetitionName = competitionName,
                                    corrected, observed, Diff = corrected-observed)]
      if(indispecies == "Jack pine" & competitionName == "Intraspecific competition"){
        output <- speciesdata
      } else {
        output <- rbind(output, speciesdata)
      }
    }
  }
  write.csv(output, file.path(workPath, "data", "areaCorrectionVerify.csv"), 
            row.names = FALSE)
} else {
  output <- read.csv(file.path(workPath, "data", "areaCorrectionVerify.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% data.table
}


output[, ':='(Species = factor(Species, levels = c("Jack pine", "Trembling aspen", "Black spruce")),
              CompetitionName = factor(CompetitionName, levels = c("Intraspecific competition", 
                                                                   "Interspecific competition")))]
for(indispecies in c("Jack pine", "Trembling aspen", "Black spruce")){
  for(competitionName in c("Intraspecific competition", "Interspecific competition")){
    thedata <- output1[Species == indispecies & CompetitionName == competitionName,]
    theana <- data.table(summary(lm(Diff~areaRatio, data = thedata))$coefficients, keep.rownames = TRUE)
    ablinetable <- data.table(Species = indispecies, CompetitionName = competitionName,
                              intercept = theana$Estimate[1], slope = theana$Estimate[2], 
                              P = theana$`Pr(>|t|)`[2])
    if(indispecies == "Jack pine" & competitionName == "Intraspecific competition"){
      ablineoutput <- ablinetable
    } else {
      ablineoutput <- rbind(ablineoutput, ablinetable)
    }
  }
}
ablineoutput[,':='(intercept = intercept/(10^floor(log10(abs(intercept)))),
                   P = P/(10^(floor(log10(abs(P)))+1)))]
ablineoutput[,':='(functionText = paste("y==", round(intercept, 2), "~+~", round(slope, 2), "~x~(italic(P)==", round(P, 2), ")", sep = ""),
                   areaRatio = 0.4, Diff = 5000)]
ablineoutput[,':='(Species = factor(Species, levels = c("Jack pine", "Trembling aspen", "Black spruce")),
                   CompetitionName = factor(CompetitionName, levels = c("Intraspecific competition", 
                                                                        "Interspecific competition")))]
titleb <- ablineoutput[CompetitionName == "Intraspecific competition" & Species == "Jack pine",][
  ,':='(areaRatio = -Inf, Diff = Inf, text = "b")]
alphaBetaTable <- data.table::copy(ablineoutput)[,.(Species, CompetitionName, areaRatio, 
                                                    Diff = -5000)]
alphaBetaTable[Species == "Jack pine" & CompetitionName == "Intraspecific competition",
               text := paste("(alpha==0.4~~beta==1.1)")]
alphaBetaTable[Species == "Jack pine" & CompetitionName == "Interspecific competition",
               text := paste("(alpha==5.3~~beta==1.0)")]
alphaBetaTable[Species == "Black spruce" & CompetitionName == "Intraspecific competition",
               text := paste("(alpha==1.0~~beta==1.0)")]
alphaBetaTable[Species == "Black spruce" & CompetitionName == "Interspecific competition",
               text := paste("(alpha==0.7~~beta==0.6)")]
alphaBetaTable[Species == "Trembling aspen" & CompetitionName == "Intraspecific competition",
               text := paste("(alpha==0.8~~beta==0.5)")]
alphaBetaTable[Species == "Trembling aspen" & CompetitionName == "Interspecific competition",
               text := paste("(alpha==4.0~~beta==0.7)")]

thefigure <- ggplot(data = output, aes(x = areaRatio, y = Diff))+
  facet_grid(Species~CompetitionName)+
  geom_point(col = "gray")+
  geom_smooth(method = "lm", se = FALSE, col = "red",
              size = 1)+
  geom_text(data = ablineoutput, aes(x = areaRatio, y = Diff, label = functionText), col = "black",
            parse = TRUE, hjust = 0)+
  geom_text(data = titleb, aes(x = areaRatio, y = Diff, label = text), vjust = 1.5, hjust = -1.5,
            size = 7)+
  geom_text(data = alphaBetaTable, aes(x = areaRatio, y = Diff, label = text), parse = TRUE,
            hjust = 0)+
  scale_y_continuous(name = "Difference between observed and corrected competition index in gray area")+
  scale_x_continuous(name = expression(paste("Ratio of ", area[ABC], " to ", area[gray])))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        strip.text = element_text(size = 15),
        strip.background = element_rect(colour = "white"))
library(png);library(grid); library(gridExtra)
img <- readPNG(file.path(workPath, "data", "areaCorrectionVerify.png"))
g <- rasterGrob(img)
g$width <- unit(1, "npc")
g$height <- unit(1, "npc")


figureA <- ggplot(data = data.frame(x = c(0, 1), y = c(0, 1)),
                  aes(x = x, y = y))+
  annotation_custom(g, xmin = 0, xmax = 1, ymin = 0, ymax = 1)+
  geom_text(data = data.frame(x = -Inf, y = Inf, text = "a"), aes(x = x, y = y, label = text),
            size = 7, vjust = 1.5, hjust = -1.5)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_blank(),
        plot.margin = margin(0,0,0,0),
        panel.margin = margin(0,0,0,0))

layout <- rbind(c(NA, 1, 1, NA, NA), rep(2, 5), rep(2, 5), rep(2, 5))
allfigure <- grid.arrange(figureA, thefigure, layout_matrix = layout)
ggsave(file.path(workPath, "TablesFigures", "areaCorrectionVerify.png"), 
       allfigure, width = 7, height = 10)





