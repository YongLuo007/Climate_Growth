rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "ClimateModelSelection168Plots.RData"))
analysesData[,SA := IniFA+2.5]


for(i in 1:length(allFixedCoeff)){
  modelName <- unlist(strsplit(names(allFixedCoeff)[i], "_", fixed = T))
  coeffIndi <- cbind(data.table(Model = rep(names(allFixedCoeff)[i], nrow(allFixedCoeff[[i]])),
                                Species = modelName[1],
                                Climate = modelName[2]), 
                     
                     allFixedCoeff[[i]])
  if(i == 1){
    allCoeff <- coeffIndi
  } else {
    allCoeff <- rbind(allCoeff, coeffIndi)
  }
}
for(i in 1:nrow(allCoeff)){
  allCoeff$rn[i] <- paste(sort(unlist(strsplit(allCoeff$rn[i], ":", fixed = TRUE))), collapse = ":")
}
OverallResults <- data.table::copy(allCoeff)
OverallResults <- OverallResults[Species != "All species",]
temperature <- c("ATA", "GSTA", "NONGSTA")
CMI <- c("ACMIA", "GSCMIA")
CO2 <- c("ACO2A")
longcol <- c(temperature, CMI, CO2)
fourspecies <- c("Jack pine", "Trembling aspen", "Black spruce")
m <- 1
for(indispecies in fourspecies){
  
  speciesData <- analysesData[Species == indispecies,]
  
  climateWithClimateTable <- rbind(data.table(expand.grid(Species = indispecies,
                                                          Climate = as.character(longcol), 
                                                          OntogenyName = "DBH",
                                                          Ontogeny = c(min(speciesData$IniDBH),
                                                                          exp(mean(log(speciesData$IniDBH))),
                                                                          max(speciesData$IniDBH)),
                                                          stringsAsFactors = FALSE)),
                                   data.table(expand.grid(Species = indispecies,
                                                          Climate = as.character(longcol), 
                                                          OntogenyName = "SA",
                                                          Ontogeny = c(min(speciesData$SA),
                                                                          exp(mean(log(speciesData$SA))),
                                                                          max(speciesData$SA)),
                                                          stringsAsFactors = FALSE))) %>% data.table
  climateWithClimateTable[OntogenyName == "DBH",
                          ':='(Ontogenyctd = log(Ontogeny)-mean(log(speciesData$IniDBH)))]
  climateWithClimateTable[OntogenyName == "SA",
                          ':='(Ontogenyctd = log(Ontogeny)-mean(log(speciesData$SA)))]
  climateWithClimateTable[,xscale:=0]
  
  climateWithClimateTable[OntogenyName == "DBH" & Ontogeny == min(speciesData$IniDBH), 
                          xscale:=as.numeric(m)]
  climateWithClimateTable[OntogenyName == "DBH" & Ontogeny == exp(mean(log(speciesData$IniDBH))), 
                          xscale:=as.numeric(m+1)]
  climateWithClimateTable[OntogenyName == "DBH" & Ontogeny == max(speciesData$IniDBH),
                          xscale:=as.numeric(m+2)]
  climateWithClimateTable[OntogenyName == "SA" & Ontogeny == min(speciesData$SA), 
                          xscale:=as.numeric(m)]
  climateWithClimateTable[OntogenyName == "SA" & Ontogeny == exp(mean(log(speciesData$SA))),
                          xscale:=as.numeric(m+1)]
  climateWithClimateTable[OntogenyName == "SA" & Ontogeny == max(speciesData$SA),
                          xscale:=as.numeric(m+2)]
  m <- m+3
  
  if(indispecies == "Jack pine"){
    alloutput <- climateWithClimateTable
  } else {
    alloutput <- rbind(alloutput, climateWithClimateTable)
  }
}


maineffectTable <- OverallResults[rn == "Climatectd",][, .(Species, Climate, mainEffect = Value,
                                                           mainEffect_SE = Std.Error)]

interactionTable <- OverallResults[rn %in% c("Climatectd:logDBHctd", "Climatectd:logSActd"),][
  ,.(rn, Species, Climate, interactEff = Value, 
     interactEff_SE = Std.Error)]
interactionTable[rn == "Climatectd:logDBHctd", OntogenyName:="DBH"]
interactionTable[rn == "Climatectd:logSActd", OntogenyName:="SA"]
interactionTable[,rn:=NULL]

alloutput1 <- dplyr::left_join(alloutput, maineffectTable, by = c("Species", "Climate")) %>%
  data.table
alloutput1[is.na(mainEffect), ':='(mainEffect = 0, mainEffect_SE = 0)]

alloutput1 <- dplyr::left_join(alloutput1, interactionTable, by = c("Species", "Climate", "OntogenyName")) %>%
  data.table
alloutput1 <- alloutput1[is.na(interactEff),':='(interactEff = 0, 
                                                 interactEff_SE = 0)]
alloutput1[,':='(Effect = mainEffect+Ontogenyctd*interactEff, 
                 Effect_SE = sqrt(mainEffect_SE^2+interactEff_SE^2))]
alloutput1[xscale %in% c(2, 5, 8, 11), ':='(Effect = mainEffect, 
                                                 Effect_SE = mainEffect_SE)]
climateWithCompTable <- data.table::copy(alloutput1)
climateWithCompTable[Climate %in% c("ATA", "ACMIA", "ACO2A") & OntogenyName == "DBH",
                     SeasonComp:="WholeDBH"]
climateWithCompTable[Climate %in% c("ATA", "ACMIA", "ACO2A") & OntogenyName == "SA",
                     SeasonComp:="WholeSA"]

climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A") & OntogenyName == "DBH",
                     SeasonComp:="GSDBH"]
climateWithCompTable[Climate %in% c("GSTA", "GSCMIA", "GSCO2A") & OntogenyName == "SA",
                     SeasonComp:="GSSA"]

climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A") & OntogenyName == "DBH", 
                     SeasonComp:="NGSDBH"]
climateWithCompTable[Climate %in% c("NONGSTA", "NONGSCMIA", "NONGSCO2A") & OntogenyName == "SA", 
                     SeasonComp:="NGSSA"]

climateWithCompTable[, SeasonComp:=factor(SeasonComp, 
                                          levels = c("WholeDBH", "WholeSA", 
                                                     "GSDBH", "GSSA", 
                                                     "NGSDBH", "NGSSA"))]

climateWithCompTable[,':='(Species = factor(Species, levels = fourspecies))]

climateWithCompTable[Climate %in% c("ATA", "GSTA", "NONGSTA"), ClimateName:="Temperature"]
climateWithCompTable[Climate %in% c("ACMIA", "GSCMIA", "NONGSCMIA"), ClimateName:="CMI"]
climateWithCompTable[Climate %in% c("ACO2A", "GSCO2A", "NONGSCO2A"), ClimateName:="CO2"]
climateWithCompTable[, ClimateName:=factor(ClimateName, levels = c("Temperature", "CMI", "CO2"))]



startPoints <- climateWithCompTable[xscale %in% c(1, 4, 7, 10),][, ':='(x = xscale+1)]
endPoints <- climateWithCompTable[xscale %in% c(3, 6, 9, 12), ][, .(Species, Climate, OntogenyName, 
                                                                    xend = xscale-1, 
                                                                    Effectend = Effect, Effectend_SE = Effect_SE)]
segmentPoints <- setkey(startPoints, Species, Climate, OntogenyName)[setkey(endPoints, Species, Climate, OntogenyName), 
                                                                        nomatch = 0]
segmentPoints[,':='(Species = factor(Species, levels = c("Jack pine", "Trembling aspen",
                                                         "Black spruce")))]
mainEffect <- climateWithCompTable[xscale %in% c(2, 5, 8, 11) & mainEffect != 0,][,':='(x = xscale)]
mainEffect[,':='(Species = factor(Species, levels = c("Jack pine", "Trembling aspen",
                                                      "Black spruce")))]


newlabels1 <- list("Temperature" = "Temperature effect",
                   "CMI" = "CMI effect",
                   'CO2'= expression(paste("C", O[2], " effect")))
newlabels2 <- list("WholeDBH" = "DBH",
                   "WholeSA" = "SA",
                   "GSDBH" = "DBH",
                   "GSSA" = "SA",
                   "NGSDBH" = "DBH",
                   "NGSSA" = "SA")

figure_labeller <- function(variable,value){
  if(variable == "ClimateName"){
    return(newlabels1[value])
  } else if (variable == "SeasonComp"){
    return(newlabels2[value])
  } 
}
majorYpanelbreaklines <- data.table(SeasonComp = factor(c("WholeDBH", "GSDBH", "NGSDBH"),
                                                        levels = c("WholeDBH", "WholeSA", 
                                                                   "GSDBH", "GSSA", 
                                                                   "NGSDBH", "NGSSA")),
                                    x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

majorXpanelbreaklines <- data.table(ClimateName = factor(c("CO2"),
                                                         levels = c("Temperature", "CMI", 
                                                                    "CO2")),
                                    x = Inf, xend = -Inf, y = -Inf, yend = -Inf)
a <- unique(segmentPoints[,.(SeasonComp, ClimateName)],
            by = c("SeasonComp", "ClimateName"))
a[,':='(x = -1, xend = Inf, Effect = 0, yend = 0)]
labeltexts <- data.table(ClimateName = factor("Temperature",
                                              levels = c("Temperature", "CMI", 
                                                         "CO2")), 
                         SeasonComp = factor(c("WholeDBH", 
                                               "GSDBH", "NGSDBH"),
                                             levels = c("WholeDBH", "WholeSA", 
                                                        "GSDBH", "GSSA", 
                                                        "NGSDBH", "NGSSA")),
                         labels = letters[1:3],
                         x = 0.5, y = Inf)

FigureB <- ggplot(data = segmentPoints[Effect != Effectend, ], aes(x = x, y = Effect))+
  facet_grid(ClimateName~SeasonComp,
             scales = "free_y", switch = "both", 
             labeller = figure_labeller, drop = FALSE)+
  geom_point(data = mainEffect, aes(x = x, y = mainEffect, group = Species, col = Species),
             size = 0.5)+
  geom_segment(data = a, aes(x = x, xend = xend, 
                             y = Effect, yend = yend), linetype = 2, col = "gray", size = 1)+
  geom_errorbar(aes(group = Species, ymin = Effect-1.98*Effect_SE, ymax = Effect+1.98*Effect_SE), 
                width = 2, col = "gray", size = 1)+
  geom_errorbar(aes(group = Species, ymin = Effectend-1.98*Effectend_SE, ymax = Effectend+1.98*Effectend_SE), 
                width = 2, col = "gray", size = 1)+
  geom_errorbar(data = mainEffect[interactEff == 0, ], aes(x = x, ymin = mainEffect-mainEffect_SE, 
                                                           ymax = mainEffect+mainEffect_SE, group = Species),
                col = "gray", width = 2, size = 1)+
  geom_segment(aes(group = Species, col = Species, xend = xend, yend = Effectend), 
               arrow = arrow(length = unit(0.1, "npc")), size = 1)+
  geom_segment(data = majorYpanelbreaklines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1.5, col = "black")+
  geom_segment(data = majorXpanelbreaklines, aes(x = x, xend = xend, y = y, yend = yend), 
               size = 1, col = "black")+
  geom_point(data = mainEffect[interactEff == 0, ], aes(x = x, y = mainEffect, group = Species, col = Species),
             show.legend = FALSE)+
  geom_text(data = labeltexts, aes(x = x, y = y, label = labels),
            vjust = 1, size = 7)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        # axis.title.x = element_text(colour = "black", size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85, 0.45),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.background = element_rect(colour = "black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 15, vjust = -2))
ggsave(file.path(workPath, "TablesFigures",
                 "Figure S5. Climate associations with ontogeny.png"), FigureB,
       width = 10, height = 7)
