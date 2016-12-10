rm(list = ls())
library(relaimpo);library(data.table);library(ggplot2); library(dplyr)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "finalYearModels.RData"))
firstRun <- FALSE
if(firstRun){
  bestIDVs <- lapply(bestFomula, function(s) as.character(s)[[3]])
  bestIDVs <- lapply(bestIDVs, function(x) {unlist(strsplit(x, " + ", fixed = TRUE))})
  bestIDVsnew <- list()
  for(i in 1:length(bestIDVs)){
    indibestIDV <- bestIDVs[[i]]
    indibestIDV <- unlist(lapply(lapply(lapply(lapply(indibestIDV, function(x) strsplit(x, ":")), function(s) unlist(s)),
                                        function(f) sort(f)), function(d) paste(d, collapse = "")))
    bestIDVsnew[[i]] <- indibestIDV
  }
  bestIDVsnew <- lapply(bestIDVsnew, function(x) gsub("ctd", "", x))
  names(bestIDVsnew) <- names(bestFomula)
  formulas <- lapply(bestIDVsnew, function(x) paste("logY~", paste(x, collapse = "+")))
  names(formulas) <- names(bestFomula)
  ontogeny <- c("logDBH", "logSA", "logDBHlogSA")
  competition <- c("logIntraH", "logInterH", "logInterHlogIntraH")
  climate <- "Year"
  ontogenyCompetition <- c("logDBHlogInterH", "logIntraHlogSA", 
                           "logInterHlogSA", "logDBHlogIntraH")
  ontogenyClimate <- c("logDBHYear","logSAYear")
  competitionClimate <- c("logIntraHYear", 
                          "logInterHYear")
  
  for(indispecies in studySpecies){
    speciesData <- analysesData[DataType == indispecies, ][
      ,':='(logY = log(BiomassGR),
            logDBH = log(IniDBH)-mean(log(IniDBH)),
            logSA = log(SA) - mean(log(SA)),
            logIntraH = log(IntraH+1) - mean(log(IntraH+1)),
            logInterH = log(InterH+1)- mean(log(InterH+1)))]
    speciesData[,':='(logDBHlogInterH = logDBH*logInterH, 
                      logDBHlogIntraH = logDBH*logIntraH,
                      logDBHlogSA = logDBH*logSA, 
                      logDBHYear = logDBH*Year, 
                      logInterHlogIntraH = logInterH*logIntraH, 
                      logInterHlogSA = logInterH*logSA,
                      logInterHYear = logInterH*Year,
                      logIntraHlogSA = logIntraH*logSA, 
                      logIntraHYear = logIntraH*Year,
                      logSAYear = logSA*Year)]
    allIDVS <- data.table(Variable = bestIDVsnew[[indispecies]])
    allIDVS[Variable %in% ontogeny, Group:="Ontogeny"]
    allIDVS[Variable %in% competition, Group:="Competition"]
    allIDVS[Variable %in% climate, Group:="Climate"]
    allIDVS[Variable %in% ontogenyCompetition, Group:="OntogenyCompetition"]
    allIDVS[Variable %in% ontogenyClimate, Group:="OntogenyClimate"]
    allIDVS[Variable %in% competitionClimate, Group:="CompetitionClimate"]
    groupnames <- unique(allIDVS$Group)
    grouplist <- list()
    for(i in groupnames){
      grouplist[[i]] <- allIDVS[Group == i,]$Variable
    }
    suppressWarnings(des <- svydesign(id=~PlotID,data = speciesData))
    bt1 <- boot.relimp(as.formula(formulas[[indispecies]]),
                       data = speciesData,
                       type = c("lmg"),
                       b = 1000,
                       design = des,
                       groups = grouplist,
                       groupnames = groupnames,
                       rela = TRUE)
    bt <- booteval.relimp(bt1, level = c(0.95))
    tempfixedLmg <- attributes(bt)
    temptable <- data.table(Species = indispecies, Group = tempfixedLmg$namen[-1],
                            importance = tempfixedLmg$lmg, 
                            importanceLower = as.numeric(tempfixedLmg$lmg.lower),
                            importanceUpper = as.numeric(tempfixedLmg$lmg.upper))
    temptable[Group %in% ontogeny, Group:="Ontogeny"]
    temptable[Group %in% competition, Group:="Competition"]
    temptable[Group %in% climate, Group:="Climate"]
    temptable[Group %in% ontogenyCompetition, Group:="OntogenyCompetition"]
    temptable[Group %in% ontogenyClimate, Group:="OntogenyClimate"]
    temptable[Group %in% competitionClimate, Group:="CompetitionClimate"]
    cat("Species: ", indispecies, "is done. \n")
    if(indispecies == "All species"){
      importanceTable <- temptable
    } else {
      importanceTable <- rbind(importanceTable, temptable)
    }
    cat("Species", indispecies, "is done. \n")
  }
  
  newimportanceTable <- importanceTable[,.(importance=importance*100, 
                                           importanceLower=importanceLower*100,
                                           importanceUpper=importanceUpper*100), 
                                        by = c("Species", "Group")]
  
  write.csv(newimportanceTable, file.path(workPath, "Results", "importanceTable.csv"), row.names = F)
} else {
  newimportanceTable <- read.csv(file.path(workPath, "Results", "importanceTable.csv"), header = TRUE,
                                 stringsAsFactors = FALSE) %>% data.table
}
newimportanceTable <- newimportanceTable[Species != "All species",]

newimportanceTable[,':='(Species = factor(Species,
                                          levels = c("Jack pine", "Other species",
                                                     "Trembling aspen",
                                                     "Black spruce")),
                         Group = factor(Group,levels = c("Ontogeny", "Competition",
                                                         "Ontogeny+Competition", "Climate", 
                                                         "Competition+Climate", "Ontogeny+Climate"),
                                        labels = c("Ontogeny", "Competition",
                                                   "Ontogeny×Competition", "Year", 
                                                   "Year×Competition",
                                                   "Year×Ontogeny")))]
subtitles <- data.table(x = 6, y = -2, 
                        Species = factor(c("Jack pine",
                                           "Trembling aspen", 
                                           "Black spruce"),
                                         levels = c("Jack pine", "Other species",
                                                    "Trembling aspen",
                                                    "Black spruce")),
                        labels = c("a", "b", "c"))

segmentstable <- rbind(c(-Inf, Inf, -Inf, -Inf),
                       c(-Inf, -Inf, -Inf, Inf))
segmentstable <- data.frame(rbind(segmentstable, segmentstable,
                                  segmentstable))
names(segmentstable) <- c("x", "xend", "y", "yend")
segmentstable$Species <- factor(sort(rep(c("Jack pine", 
                                           "Trembling aspen",
                                           "Black spruce"), 2)),
                                levels=c("Jack pine", "Other species",
                                         "Trembling aspen", "Black spruce"))
for(i in 2:4){
  label1 <- data.table(Species = names(allFixedCoeff[i]), 
                       marR2 = paste("Marginal~R^2==", round(allFixedCoeff[[i]]$marR2[1], 2)),
                       conR2 = paste("Conditional~R^2==", round(allFixedCoeff[[i]]$conR2[1], 2)))
  if(i == 2){
    labelAll <- label1
  } else {
    labelAll <- rbind(labelAll, label1)
  }
}
labelAll[, ':='(y = 8, xmar = 3.6, xcon = 3)]
labelAll[, Species := factor(Species, levels = c("Jack pine", "Other species",
                                                 "Trembling aspen",
                                                 "Black spruce"))]
importanceFigure <- ggplot(data = newimportanceTable[Species != "Other species"],
                           aes(x = Group, y = importance))+
  geom_bar(aes(col = Group, fill = Group), stat = 'identity')+
  geom_errorbar(aes(ymin = importanceLower, ymax = importanceUpper), col = "gray", width = 0.25)+
  scale_x_discrete(limits = rev(c("Ontogeny", "Competition",
                                  "Ontogeny×Competition", "Year", 
                                  "Year×Competition",
                                  "Year×Ontogeny")))+
  scale_y_continuous(name = "Relative importance (%)")+
  
  geom_segment(data = segmentstable, 
               aes(x = x, xend = xend, y = y, yend = yend))+
  geom_text(data = subtitles, aes(x = x, y = y, label = labels), size = 10)+
  geom_text(data = labelAll, aes(x = xmar, y = y, label = marR2), 
            size = 5, parse = TRUE, hjust = 0)+
  geom_text(data = labelAll, aes(x = xcon, y = y, label = conR2), 
            size = 5, parse = TRUE, hjust = 0)+
  facet_wrap(~Species, ncol = 2, drop = FALSE)+
  scale_fill_manual(name = "Predictor group", values = c("cyan", "darkgreen", "darkblue", 
                                                         "purple", "goldenrod", "brown"),
                    labels = c("Ontogeny", "Competition",
                               "Ontogeny\n×Competition", "Year", 
                               "Year\n×Competition",
                               "Year\n×Ontogeny"))+
  scale_color_manual(guide = "none", values = c("cyan", "darkgreen", "darkblue", 
                                                "purple", "goldenrod", "brown"))+
  guides(guide_legend(nrow = 3))+ 
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank(),
        legend.key.height = unit(1.5, "lines"),
        legend.title = element_text(13),
        legend.text = element_text(size = 10),
        legend.position = c(0.8, 0.8),
        legend.background = element_rect(colour = "black"),
        strip.text = element_blank())

ggsave(file = file.path(workPath, "TablesFigures", "Figure 2. RelativeImportance.png"),
       importanceFigure,
       width = 9, height = 7.5)

