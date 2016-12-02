rm(list = ls())
library(SpaDES); library(ggplot2); library(grid); library(gridExtra)
library(dplyr);library(data.table)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "newAllDataRescaledComp.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table

normalityCheckData <- data.table::copy(analysesData)
normalityCheckData[,DataType:=factor(DataType, levels = c("All species", "Jack pine", "", "Trembling aspen", 
                                               "Black spruce", "Other species"))]


# visual check
normalityCheckFigure <- ggplot(normalityCheckData, aes(x = BiomassGR))+
  geom_histogram(aes(log(BiomassGR)), bins = 40)+
  facet_wrap(~DataType, ncol = 3, scales = "free_y", drop = FALSE)+
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
studyspecies <- c("All species", "Jack pine", "Trembling aspen", "Black spruce", "Other species")
for(indispecies in studyspecies){
  independentV1 <- data.frame(analysesData[DataType == indispecies, .(SA,
                                                                      Year)])
  
  corTable <- data.table(Species = indispecies,cor(independentV1), keep.rownames = T, key = "rn")
  viftest <- data.table(corvif(independentV1), keep.rownames = TRUE, key = "rn")
  viftest <- corTable[viftest, nomatch = 0]
  setnames(viftest, "rn", "Variable")
  if(indispecies == "All species"){
    viftestoutput <- viftest
  } else {
    viftestoutput <- rbind(viftestoutput, viftest)
  }
}
viftestoutput <- cbind(viftestoutput[,.(Species, Variable)], 
                       round(viftestoutput[,.(DBH, SA, IntraH, InterH, Year, VIF=GVIF)], 2))

print(viftestoutput)
#     Species Variable   DBH InterH IntraH   RBI  Year  VIF
# 1:     All      DBH  1.00  -0.32  -0.51  0.57  0.00 1.73
# 2:     All   InterH -0.32   1.00   0.22 -0.35 -0.02 1.17
# 3:     All   IntraH -0.51   0.22   1.00 -0.43  0.06 1.42
# 4:     All      RBI  0.57  -0.35  -0.43  1.00 -0.07 1.61
# 5:     All     Year  0.00  -0.02   0.06 -0.07  1.00 1.01
# 6:      JP      DBH  1.00  -0.14  -0.58  0.47 -0.12 1.66
# 7:      JP   InterH -0.14   1.00   0.05 -0.02  0.04 1.02
# 8:      JP   IntraH -0.58   0.05   1.00 -0.48  0.22 1.72
# 9:      JP      RBI  0.47  -0.02  -0.48  1.00 -0.07 1.40
# 10:      JP     Year -0.12   0.04   0.22 -0.07  1.00 1.05
# 11:      TA      DBH  1.00  -0.28  -0.66  0.49  0.06 2.28
# 12:      TA   InterH -0.28   1.00   0.01 -0.12  0.13 1.19
# 13:      TA   IntraH -0.66   0.01   1.00 -0.45  0.01 1.99
# 14:      TA      RBI  0.49  -0.12  -0.45  1.00 -0.05 1.38
# 15:      TA     Year  0.06   0.13   0.01 -0.05  1.00 1.05
# 16:      BS      DBH  1.00  -0.39  -0.50  0.63 -0.01 1.91
# 17:      BS   InterH -0.39   1.00   0.24 -0.40 -0.06 1.25
# 18:      BS   IntraH -0.50   0.24   1.00 -0.45  0.03 1.39
# 19:      BS      RBI  0.63  -0.40  -0.45  1.00 -0.12 1.86
# 20:      BS     Year -0.01  -0.06   0.03 -0.12  1.00 1.03
# 21:   Other      DBH  1.00  -0.44  -0.40  0.66  0.08 2.09
# 22:   Other   InterH -0.44   1.00   0.33 -0.39  0.11 1.35
# 23:   Other   IntraH -0.40   0.33   1.00 -0.34 -0.14 1.29
# 24:   Other      RBI  0.66  -0.39  -0.34  1.00 -0.13 1.95
# 25:   Other     Year  0.08   0.11  -0.14 -0.13  1.00 1.13

# the vif is not an issue for analyses, as indicated by the GVIF value smaller than 4


##
### check the linearity between dependent variable and independent variable visually
###
analysesDataLongForm <- reshape(data = analysesData[,.(Species, BiomassGR, Year, 
                                                       logIntraH = log(Hegyi*IntraHegyiRatio+1),
                                                       logInterH = log(Hegyi*(1-IntraHegyiRatio)+1), 
                                                       RBI, logDBH = log(IniDBH))], 
                                varying = c("Year", "logIntraH", "logInterH", "RBI", "logDBH"),
                                v.names = "Value",
                                timevar = "IndependentV",
                                times = c("Year", "logIntraH", "logInterH", "RBI", "logDBH"),
                                direction = "long")
analysesDataLongForm[!(Species %in% majorspecies), Species:="Other"]
analysesDataLongFormall <- data.table::copy(analysesDataLongForm)[,Species:="All"]
analysesDataLongForm <- rbindlist(list(analysesDataLongForm, analysesDataLongFormall))
rm(analysesDataLongFormall)
analysesDataLongForm[,Species:=factor(Species, levels = c("All", "JP", "TA", "BS", "Other"),
                                      labels = c("All species", "Jack pine", "Trembling aspen",
                                                 "Black spruce", "Other species"))]

figure_linearity <- ggplot(data = analysesDataLongForm, aes(y = log(BiomassGR)))+
  geom_point(aes(x = Value, group = Species))+
  facet_grid(Species~IndependentV, scales = "free")+
  geom_segment(aes(x = -Inf, xend = Inf, y = -Inf, yend = -Inf), size = 1)+
  geom_segment(aes(x = -Inf, xend = -Inf, y = -Inf, yend = Inf), size = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))


dev(4)
clearPlot()
ggsave(file = file.path(workPath, "TablesFigures", "Figure S3 LinearityCheck.png"), 
       figure_linearity, width = 14, height = 9)





