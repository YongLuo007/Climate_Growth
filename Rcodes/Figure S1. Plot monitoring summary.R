rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(gridExtra)

workPath <- "~/GitHub/Climate_Growth"

analysesData <- fread(file.path(workPath, "data", "Year10Analyses",
                                   "finalData10.csv"))
analysesData <- analysesData[allCensusLiveTree == "yes",]

allspeciesdata <- data.table::copy(analysesData)
FigureS1Data <- data.table(PlotID = character(), Year = numeric(), IniYear = numeric(),
                           FinYear = numeric())
plots <- unique(allspeciesdata$PlotID)
# measurementLength summary
measurementLengthData <- unique(allspeciesdata[,.(PlotID, IniYear, FinYear)],
                                by = c("PlotID", "IniYear"))[,':='(Length = FinYear-IniYear,
                                                                   Year = (FinYear+IniYear)/2)]
measurementLengthData[,':='(FirstYear = min(IniYear), LastYear = max(FinYear)), by = PlotID]
newplotiddata <- unique(measurementLengthData[,.(PlotID, FirstYear, LastYear)], by = "PlotID")
newplotiddata <- newplotiddata[order(FirstYear, LastYear),]
newplotiddata[, newPlotID:=1:nrow(newplotiddata)]
measurementLengthData <- setkey(measurementLengthData, PlotID)[setkey(newplotiddata[,.(PlotID, newPlotID)], PlotID), nomatch = 0]
Non10YearsData <- measurementLengthData[Length != 10,]

subFigureS1 <- reshape(data = as.data.frame(newplotiddata), varying = c("FirstYear", "LastYear"), 
                       v.names = "Year",
                       direction = "long", times = c("FirstYear", "LastYear"),
                       timevar = "FirstOrLast") %>%
  data.table
subFigureS1_1 <- subFigureS1[,.(Year = mean(Year), sdYear = sd(Year),
                                y = 130), by = FirstOrLast]

FigureS1_a <- ggplot(data = subFigureS1, aes(x = Year, y = newPlotID))+
  geom_line(aes(group = newPlotID), size = 0.5, col = "gray")+
  geom_point(data = subFigureS1_1, aes(x = Year, y = y), size = 2)+
  geom_errorbarh(data = subFigureS1_1, aes(y = y, xmin = Year-sdYear, xmax = Year+sdYear),
                 size = 0.5, height = 7)+
  geom_segment(data = Non10YearsData, aes(y = newPlotID, x = IniYear, xend = FinYear, yend = newPlotID,
                                         group = newPlotID, col = as.factor(Length)))+
  scale_color_manual(name = "Measurement length", values = c("red", "green"),
                     label = paste(c(8, 10), " years"))+
  annotate("text", x = subFigureS1_1$Year, y = 138, size = 4,
           label = c("Year at first census (mean±SD)", 
                     "Year at last census (mean±SD)"))+
  scale_y_continuous(name = "Plot identity", limits = c(0, 140), breaks = seq(0, 180, by = 30))+
  scale_x_continuous(name = "Year", limits = c(1985, 2012), breaks = seq(1985, 2011, by = 5))+
  annotate("text", x = 1985, y = 140, label = "a", size = 11)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = c(0.15, 0.65),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))


allspeciesdatanew <- data.table::copy(analysesData)


NofTreeData <- allspeciesdatanew[,.(NofTree=length(unique(uniTreeID)),
                                 IniYear = min(IniYear),
                                 FinYear = max(FinYear)), by = c("PlotID", "Species_Group")]
figureS1_bData <- data.table(Species = character(), NofTree = numeric(), Year = numeric())
for(i in 1:nrow(NofTreeData)){
  adddata <- data.table(Species = NofTreeData$Species_Group[i], 
                        NofTree = NofTreeData$NofTree[i],
                        Year = seq(NofTreeData$IniYear[i], NofTreeData$FinYear[i], by = 1))
  figureS1_bData <- rbind(figureS1_bData, adddata)
}
figureS1_bData <- figureS1_bData[,.(NofTree = sum(NofTree)), by = c("Year", "Species")]
figureS1_bData[Species == "Minor species", Species:="Minor species group"]
studySpecies <- c("All species", "Jack pine", "Trembling aspen",
                  "Black spruce", "Minor species group")
figureS1_bData[,Species:=factor(Species, 
                                levels = studySpecies)]
summaryPlot <- allspeciesdatanew[,.(NofPlot = length(unique(PlotID)), NofTree = length(unique(uniTreeID)),
                                    NofObs = length(IniBA)), by = Species_Group]
summaryPlot[, Species := Species_Group]
summaryPlot[Species_Group == "Minor species", Species:="Minor species group"]
summaryPlot[Species_Group == "All species", Species:="All trees"]
# summaryPlot[, Species:=factor(Species, levels = studySpecies)]
summaryPlot[,legendTexts:=paste("\n", Species, ": \n", NofPlot, " plots, ", NofTree, " trees, ",
                                NofObs, " observations.", sep = "")]
studySpecies1 <- data.table(Species = c("All trees", "Jack pine", 
                                        "Trembling aspen",
                  "Black spruce", "Minor species group"),
                  orders = 1:5)
a <- dplyr::left_join(studySpecies1, summaryPlot, by = "Species")

FigureS1_b <- ggplot(data = figureS1_bData, aes(x = Year, y = NofTree))+
  geom_line(aes(group = Species, col = Species), size = 1)+
  scale_y_continuous(name = "Number of tree", limits = c(0, 16000), breaks = seq(0, 10000, by = 2000))+
  scale_x_continuous(name = "Year", limits = c(1985, 2012), breaks = seq(1985, 2011, by = 5))+
  scale_color_manual(values = c("white", "red", "blue", "green", "pink"),
                     labels = a$legendTexts)+
  guides(colour  = guide_legend(title = NULL, nrow = 3))+
  annotate("text", x = 1985, y = 16000, label = "b", size = 11)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = c(0.53, 0.88),
        legend.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.key = element_rect(colour = "white"),
        # legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

FigureS1_a_Grob <- ggplotGrob(FigureS1_a)
FigureS1_b_Grob <- ggplotGrob(FigureS1_b)


FigureS1_a_Grob$widths <- FigureS1_b_Grob$widths
FigureS1_a_Grob$heights <- FigureS1_b_Grob$heights
dev(4)
clearPlot()
plotlayout <- rbind(c(1), c(2))
FigureS1 <- grid.arrange(FigureS1_a_Grob, FigureS1_b_Grob, 
                  layout_matrix = plotlayout)
ggsave(file = file.path(workPath, "TablesFigures", "Figure S1. MonitoringSummary.png"), FigureS1,
       width = 9, height = 7)


