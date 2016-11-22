rm(list = ls())
library(relaimpo);library(data.table);library(ggplot2)
load(file.path("~/GitHub/Climate_Growth/Results/finalYearModels.RData"))
bestIDVsnew <- list()
for(i in 1:length(bestIDVs)){
  indibestIDV <- bestIDVs[[i]]
  indibestIDV <- unlist(lapply(lapply(lapply(lapply(indibestIDV, function(x) strsplit(x, ":")), function(s) unlist(s)),
                               function(f) sort(f)), function(d) paste(d, collapse = "")))
  bestIDVsnew[[i]] <- indibestIDV
}

lapply(bestIDVsnew, function(x) paste(x, collapse = "+"))

for(indispecies in studySpecies){
  speciesData <- finalAnalysesData[[indispecies]][
    ,':='(logDBHctdlogInterHctd = logDBHctd*logInterHctd, 
          logDBHctdlogIntraHctd = logDBHctd*logIntraHctd,
          logDBHctdlogSActd = logDBHctd*logSActd, 
          logDBHctdlogSBctd = logDBHctd*logSBctd,
          logDBHctdRBIctd = logDBHctd*RBIctd,
          logDBHctdYearctd = logDBHctd*Yearctd, 
          logInterHctdlogIntraHctd = logInterHctd*logIntraHctd, 
          logInterHctdlogSActd = logInterHctd*logSActd,
          logInterHctdlogSBctd = logInterHctd*logSBctd, 
          logInterHctdRBIctd = logInterHctd*RBIctd,
          logInterHctdYearctd = logInterHctd*Yearctd,
          logIntraHctdlogSActd = logIntraHctd*logSActd, 
          logIntraHctdlogSBctd = logIntraHctd*logSBctd,
          logIntraHctdRBIctd = logIntraHctd*RBIctd,
          logIntraHctdYearctd = logIntraHctd*Yearctd, 
          logSActdlogSBctd = logSActd*logSBctd,
          logSActdRBIctd = logSActd*RBIctd,
          logSActdYearctd = logSActd*Yearctd,
          logSBctdRBIctd = logSBctd*RBIctd, 
          logSBctdYearctd = logSBctd*Yearctd,
          RBIctdYearctd = RBIctd*Yearctd)]
  
  if(indispecies == "All species"){
    bt <- calc.relimp(logY~logDBHctd+logIntraHctd+logInterHctd+logSActd+logDBHctdYearctd+logDBHctdlogIntraHctd+logDBHctdlogInterHctd+logDBHctdlogSActd+logIntraHctdYearctd+logInterHctdYearctd+logInterHctdlogIntraHctd+logIntraHctdlogSActd+logInterHctdlogSActd,
                      data = speciesData,
                      type = c("lmg"),
                      rela = TRUE)
  } else if (indispecies == "Jack pine"){
    bt <- calc.relimp(logY~logDBHctd+logIntraHctd+logInterHctd+logSActd+logDBHctdlogIntraHctd+logDBHctdlogSActd+logIntraHctdYearctd+logSActdYearctd+logIntraHctdlogSActd+logInterHctdlogSActd,
                      data = speciesData, rela = TRUE,
                      type = c("lmg"))
  } else if (indispecies == "Trembling aspen"){
    bt <- calc.relimp(logY~logDBHctd+logIntraHctd+logInterHctd+logSActd+logDBHctdYearctd+logSActdYearctd+logInterHctdlogIntraHctd+logInterHctdlogSActd,
                      data = speciesData, rela = TRUE,
                      type = c("lmg"))
  } else if (indispecies == "Black spruce"){
    bt <- calc.relimp(logY~logDBHctd+Yearctd+logIntraHctd+logInterHctd+logSActd+logDBHctdYearctd+logDBHctdlogIntraHctd+logDBHctdlogInterHctd+logDBHctdlogSActd+logIntraHctdYearctd+logInterHctdYearctd+logInterHctdlogSActd,
                      data = speciesData, rela = TRUE,
                      type = c("lmg"))
  } else if (indispecies == "Other species"){
    bt <- calc.relimp(logY~logDBHctd+Yearctd+logIntraHctd+logInterHctd+logSActd+logDBHctdlogInterHctd+logSActdYearctd+logInterHctdlogIntraHctd+logInterHctdlogSActd,
                      data = speciesData, rela = TRUE,
                      type = c("lmg"))
  }
  tempfixedLmg <- attributes(bt)
  temptable <- data.table(Species = indispecies, Variable = tempfixedLmg$namen[-1],
                          lmg = tempfixedLmg$lmg)
  if(indispecies == "All species"){
    importanceTable <- temptable
  } else {
    importanceTable <- rbind(importanceTable, temptable)
  }
  cat("Species", indispecies, "is done. \n")
}
ontogeny <- c("logDBHctd", "logSActd", "logDBHctdlogSActd")
competition <- c("logIntraHctd", "logInterHctd", "RBIctd", "logSBctd",
                 "logInterHctdlogIntraHctd", "logIntraHctdRBIctd",
                 "logInterHctdRBIctd","logInterHctdlogSBctd",
                 "logIntraHctdlogSBctd", "logSBctdRBIctd")
climate <- "Yearctd"
ontogenyCompetition <- c("logDBHctdlogInterHctd", "logDBHctdRBIctd", 
                         "logDBHctdlogSBctd","logIntraHctdlogSActd", 
                         "logInterHctdlogSActd", "logSActdRBIctd", 
                         "logSActdlogSBctd","logDBHctdlogIntraHctd")
ontogenyClimate <- c("logDBHctdYearctd","logSActdYearctd")
competitionClimate <- c( "RBIctdYearctd", "logSBctdYearctd", 
                         "logIntraHctdYearctd", 
                         "logInterHctdYearctd")

importanceTable[Variable %in% ontogeny, Group:="Ontogeny"]
importanceTable[Variable %in% competition, Group:="Competition"]
importanceTable[Variable %in% climate, Group:="Climate"]
importanceTable[Variable %in% ontogenyCompetition, Group:="Ontogeny+Competition"]
importanceTable[Variable %in% ontogenyClimate, Group:="Ontogeny+Climate"]
importanceTable[Variable %in% competitionClimate, Group:="Competition+Climate"]


newimportanceTable <- importanceTable[,.(importance = 100*sum(lmg)), 
                                      by = c("Species", "Group")]

newimportanceTable <- rbind(newimportanceTable, 
                            data.table::copy(newimportanceTable)[1,][,Species:=" "])
newimportanceTable[,':='(Species = factor(Species,
                                          levels = c("All species", " ", "Jack pine", 
                                                     "Trembling aspen", "Black spruce", 
                                                     "Other species")),
                         Group = factor(Group,
                                        levels = c("Ontogeny", "Competition",
                                                   "Climate", "Ontogeny+Competition",
                                                   "Ontogeny+Climate", "Competition+Climate")))]
subtitles <- data.table(x = 6, y = -2, 
                        Species = factor(c("All species", " ", "Jack pine", 
                                           "Trembling aspen", "Black spruce", 
                                           "Other species"), 
                                         levels=c("All species", " ", "Jack pine", 
                                                  "Trembling aspen", "Black spruce", 
                                                  "Other species")),
                        labels = c("a", "", "b", "c", "d", "e"))

segmentstable <- rbind(c(-Inf, Inf, -Inf, -Inf),
                       c(-Inf, -Inf, -Inf, Inf))
segmentstable <- data.frame(rbind(segmentstable, segmentstable, segmentstable,
                       segmentstable, segmentstable))
names(segmentstable) <- c("x", "xend", "y", "yend")
segmentstable$Species <- factor(sort(rep(c("All species", "Jack pine", 
                                    "Trembling aspen", "Black spruce", 
                                    "Other species"), 2)),
                                levels=c("All species", " ", "Jack pine", 
                                         "Trembling aspen", "Black spruce", 
                                         "Other species"))

importanceFigure <- ggplot(data = newimportanceTable[Species != " "], aes(x = Group, y = importance))+
  geom_bar(aes(col = Group, fill = Group), stat = 'identity')+
  scale_x_discrete(limits = rev(c("Ontogeny", "Competition",
                                       "Climate", "Ontogeny+Competition",
                                       "Ontogeny+Climate", "Competition+Climate")))+
  scale_y_continuous(name = "Relative importance (%)")+
  geom_segment(data = segmentstable, 
               aes(x = x, xend = xend, y = y, yend = yend))+
  geom_text(data = subtitles, aes(x = x, y = y, label = labels), size = 10)+
  
  facet_wrap(~Species, ncol = 2)+
  coord_flip()+
  guides(col = guide_legend(title = "Variable group", nrow = 3),
         fill = guide_legend(title = "Variable group", nrow = 3))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.position = c(0.75, 0.85),
        strip.text = element_blank())

ggsave(file = file.path("~/GitHub/Climate_Growth/TablesFigures/RelativeImportance.png"),
       importanceFigure,
       width = 9, height = 7.5)

