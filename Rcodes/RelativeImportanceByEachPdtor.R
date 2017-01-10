rm(list = ls())
library(relaimpo);library(data.table);library(ggplot2); library(dplyr)
workPath <- "~/GitHub/Climate_Growth"

newimportanceTable <- read.csv(file.path(workPath, "Results", "importanceTable_indipdtor.csv"), header = TRUE,
                               stringsAsFactors = FALSE) %>% data.table
newimportanceTable[, Group:=gsub("log", "", Group)]
predictorOrder <- c("DBH", "SA", "DBH×SA", "IntraH", "InterH", "IntraH×InterH",
                    "DBH×IntraH", "DBH×InterH", "SA×IntraH", "SA×InterH", 
                    "Year", "Year×IntraH", "Year×InterH",  
                    "Year×DBH", "Year×SA")
newimportanceTable[,':='(Species = factor(Species,
                                          levels = c("Jack pine", "Other species",
                                                     "Trembling aspen",
                                                     "Black spruce")),
                         Group = factor(Group,levels = c("DBH", "SA", "DBHSA", "IntraH", "InterH", "InterHIntraH",
                                                         "DBHIntraH", "DBHInterH", "IntraHSA", "InterHSA", 
                                                         "Year", "IntraHYear", "InterHYear",  
                                                         "DBHYear", "SAYear"),
                                        labels = predictorOrder))]

newimportanceTable[Group %in% c("DBH", "SA", "DBH×SA"), newGroup := "Ontogeny"]
newimportanceTable[Group %in% c("IntraH", "InterH", "IntraH×InterH"), newGroup := "Competition"]
newimportanceTable[Group %in% c("DBH×IntraH", "DBH×InterH", "SA×IntraH", "SA×InterH"), 
                   newGroup := "Ontogeny×Competition"]
newimportanceTable[Group %in% c("Year"), newGroup := "Year"]
newimportanceTable[Group %in% c("Year×IntraH", "Year×InterH"), newGroup := "Year×Competition"]
newimportanceTable[Group %in% c("Year×DBH", "Year×SA"), newGroup := "Year×Ontogeny"]
newimportanceTable[,newGroup:=factor(newGroup, levels = c("Ontogeny", "Competition", "Ontogeny×Competition",
                                                          "Year", "Year×Competition", "Year×Ontogeny"))]


subtitles <- data.table(x = 15, y = -2, 
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

importanceFigure <- ggplot(data = newimportanceTable[Species != "Other species"],
                           aes(x = Group, y = importance))+
  facet_wrap(~Species, ncol = 2, drop = FALSE)+
  geom_bar(aes(col = newGroup, fill = newGroup), stat = 'identity')+
  geom_errorbar(aes(ymin = importanceLower, ymax = importanceUpper), col = "gray", width = 0.25)+
  scale_x_discrete(limits = rev(predictorOrder))+
  scale_y_continuous(name = "Relative importance (%)")+
  geom_segment(data = segmentstable, 
               aes(x = x, xend = xend, y = y, yend = yend))+
  geom_text(data = subtitles, aes(x = x, y = y, label = labels), size = 10)+
  scale_fill_manual(name = "Predictor group", values = c("cyan","darkgreen", "darkblue", 
                                                         "purple", "goldenrod", "brown"))+
  scale_color_manual(name = "Predictor group", values = c("cyan","darkgreen", "darkblue", 
                                "purple", "goldenrod", "brown"))+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.position = c(0.70, 0.8),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.background = element_rect(colour = "black"),
        strip.text = element_blank())

ggsave(file = file.path(workPath, "TablesFigures", "Figure S. RelativeImportanceByEachPdtor.png"),
       importanceFigure,
       width = 8, height = 10)

