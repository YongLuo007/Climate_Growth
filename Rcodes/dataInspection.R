rm(list = ls())
library(SpaDES); library(ggplot2); library(grid); library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
majorspecies <- c("JP", "TA", "BS")

normalityCheckData <- analysesData[,.(BiomassGR, Species, newspecies = Species)]
normalityCheckData[!(Species %in% majorspecies), newspecies:="Other"]
normalityCheckData <- rbind(data.table::copy(normalityCheckData)[,newspecies:="All"],
                            normalityCheckData)
normalityCheckData[,Species:=factor(newspecies, levels = c("All", "JP", "", "TA", "BS", "Other"), 
                                    labels = c("All species", "Jack pine", "", "Trembling aspen", 
                                               "Black spruce", "Other species"))]


# visual check
normalityCheckFigure <- ggplot(normalityCheckData, aes(x = BiomassGR))+
  geom_histogram(aes(log(BiomassGR)), bins = 40)+
  facet_wrap(~Species, ncol = 3, scales = "free_y", drop = FALSE)+
  scale_y_continuous(name = "Number of observations")+
  scale_x_continuous(name = "log(ABGR)")+
  geom_segment(aes(x = -Inf, xend = Inf, y = -Inf, yend = -Inf), size = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15, hjust = 0.1),
        strip.background = element_blank())

dev(4)
clearPlot()
ggsave(file = file.path(workPath, "TablesFigures", "Figure S2 NormalityCheck.png"), 
       normalityCheckFigure, width = 10, height = 10)

# check the multilinearity among independent variables using Zuur 2010 method
source(file.path(workPath, "Rcodes", "Rfunctions", "HighstatLib_Zuur.R"))
studyspecies <- c("All", "JP", "TA", "BS", "Other")
for(indispecies in studyspecies){
  if(indispecies == "All"){
    independentV1 <- data.frame(analysesData[,.(RBI, H = Hegyi, DBH=IniDBH, Year)])
  } else if (indispecies == "Other"){
    independentV1 <- data.frame(analysesData[!(Species %in% majorspecies),
                                             .(RBI, H = Hegyi, DBH=IniDBH, Year)])
  } else {
    independentV1 <- data.frame(analysesData[Species == indispecies,
                                             .(RBI, H = Hegyi, DBH=IniDBH, Year)])
  }
  corTable <- data.table(Species = indispecies,cor(independentV1), keep.rownames = T, key = "rn")
  viftest <- data.table(corvif(independentV1), keep.rownames = TRUE, key = "rn")
  viftest <- corTable[viftest, nomatch = 0]
  setnames(viftest, "rn", "Variable")
  if(indispecies == "All"){
    viftestoutput <- viftest
  } else {
    viftestoutput <- rbind(viftestoutput, viftest)
  }
}
viftestoutput <- cbind(viftestoutput[,.(Species, Variable)], 
                       round(viftestoutput[,.(DBH,H, RBI, Year, VIF=GVIF)], 2))

print(viftestoutput)
#     Species Variable   DBH     H   RBI  Year  VIF
# 1:     All      DBH  1.00 -0.54  0.57  0.00 1.69
# 2:     All        H -0.54  1.00 -0.50  0.03 1.54
# 3:     All      RBI  0.57 -0.50  1.00 -0.07 1.61
# 4:     All     Year  0.00  0.03 -0.07  1.00 1.01
# 5:      JP      DBH  1.00 -0.59  0.47 -0.12 1.65
# 6:      JP        H -0.59  1.00 -0.48  0.22 1.74
# 7:      JP      RBI  0.47 -0.48  1.00 -0.07 1.39
# 8:      JP     Year -0.12  0.22 -0.07  1.00 1.05
# 9:      TA      DBH  1.00 -0.70  0.49  0.06 2.22
# 10:      TA        H -0.70  1.00 -0.43  0.08 2.05
# 11:      TA      RBI  0.49 -0.43  1.00 -0.05 1.35
# 12:      TA     Year  0.06  0.08 -0.05  1.00 1.04
# 13:      BS      DBH  1.00 -0.57  0.63 -0.01 1.89
# 14:      BS        H -0.57  1.00 -0.54 -0.02 1.61
# 15:      BS      RBI  0.63 -0.54  1.00 -0.12 1.86
# 16:      BS     Year -0.01 -0.02 -0.12  1.00 1.03
# 17:   Other      DBH  1.00 -0.52  0.66  0.08 2.09
# 18:   Other        H -0.52  1.00 -0.45  0.02 1.40
# 19:   Other      RBI  0.66 -0.45  1.00 -0.13 1.95
# 20:   Other     Year  0.08  0.02 -0.13  1.00 1.07

# the vif is not an issue for analyses, as indicated by the GVIF value smaller than 4


##
### check the linearity between dependent variable and independent variable visually
###
analysesDataLongForm <- reshape(data = analysesData[,.(Species, BiomassGR, Year, Hegyi, RBI,
                                                       IniDBH)], 
                                varying = c("Year", "Hegyi", "RBI", "IniDBH"),
                                v.names = "Value",
                                timevar = "IndependentV",
                                times = c("Year", "Hegyi", "RBI", "IniDBH"),
                                direction = "long")


figure_linearity_1 <- ggplot(data = analysesData, aes(y = log(BiomassGR)))+
  geom_point(aes(x = Year, group = Species))+
  facet_wrap(~Species)
figure_linearity_2 <- ggplot(data = analysesData, aes(y = log(BiomassGR)))+
  geom_point(aes(x = log(Hegyi), group = Species))+
  facet_wrap(~Species)
figure_linearity_3 <- ggplot(data = analysesData, aes(y = log(BiomassGR)))+
  geom_point(aes(x = log(IniDBH), group = Species))+
  facet_wrap(~Species)
figure_linearity_4 <- ggplot(data = analysesData, aes(y = log(BiomassGR)))+
  geom_point(aes(x = RBI, group = Species))+
  facet_wrap(~Species)


figure_linearity_1G <- ggplotGrob(figure_linearity_1)
figure_linearity_2G <- ggplotGrob(figure_linearity_2)
figure_linearity_3G <- ggplotGrob(figure_linearity_3)
figure_linearity_4G <- ggplotGrob(figure_linearity_4)

figure_linearity_1G$heights <- figure_linearity_4G$heights
figure_linearity_1G$widths <- figure_linearity_4G$widths
figure_linearity_2G$heights <- figure_linearity_4G$heights
figure_linearity_2G$widths <- figure_linearity_4G$widths
figure_linearity_3G$heights <- figure_linearity_4G$heights
figure_linearity_3G$widths <- figure_linearity_4G$widths
dev(4)
clearPlot()
plotlayout <- rbind(c(1), c(2), c(3), c(4))
gridExtra::grid.arrange(figure_linearity_1G, figure_linearity_2G,
                   figure_linearity_3G, figure_linearity_4G,
                  layout_matrix = plotlayout)



