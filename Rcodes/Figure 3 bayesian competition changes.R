rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "Bayesian_YearModels.RData"))
rm(indispecies)



for(indispecies in c("All species", "Jack pine", "Trembling aspen",
                     "Black spruce", "Other species")){
  Year <- seq(min(Figure3Data[Species == indispecies, ]$Year), 
              max(Figure3Data[Species == indispecies, ]$Year), length = 50)
  intrahtrend <- data.table(Species = indispecies,
                            Competition = "IntraH",
                            Year = Year, 
                            Yearctd = Year-
                              mean(Figure3Data[Species == indispecies &Main == 1,]$Year))
  intrahtrend[,':='(Value = Yearctd*fixedEffect[Species == indispecies & Variable == "c7"]$Mean+
                      fixedEffect[Species == indispecies & Variable == "b2"]$Mean,
                    lineType = 1)]
  if(fixedEffect[Species == indispecies & Variable == "c7"]$Lower95*
     fixedEffect[Species == indispecies & Variable == "c7"]$Upper <= 0){
    intrahtrend[, lineType := 2]
  }
  
  interhtrend <- data.table(Species = indispecies,
                            Competition = "InterH",
                            Year = Year, 
                            Yearctd = Year-
                              mean(Figure3Data[Species == indispecies &Main == 1,]$Year))
  interhtrend[,':='(Value = Yearctd*fixedEffect[Species == indispecies & Variable == "c9"]$Mean+
                      fixedEffect[Species == indispecies & Variable == "b3"]$Mean,
                    lineType = 1)]
  if(fixedEffect[Species == indispecies & Variable == "c9"]$Lower95*
     fixedEffect[Species == indispecies & Variable == "c9"]$Upper <= 0){
    interhtrend[, lineType := 2]
  }
  CompetitionTrends_indis <- rbind(intrahtrend, interhtrend)
  if(indispecies == "All species"){
    CompetitionTrends <- CompetitionTrends_indis
  } else{
    CompetitionTrends <- rbind(CompetitionTrends, CompetitionTrends_indis)
  }
}
rm(intrahtrend, interhtrend, Year, indispecies)

CompetitionTrends[, ':='(Direction = factor(Competition, levels = c("IntraH", "InterH")),
                   Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen", 
                                                        "Black spruce", "Other species")),
                   lineType = factor(lineType, levels = c(1, 2)))]

equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    round(seq(min(x)+d, max(x)-d, length=n), 2)
  }
}  
yaxisLineData <- data.table(expand.grid(Direction = c("IntraH", "InterH"), 
                                        Species = c("All species", "Jack pine", "Trembling aspen", 
                                                    "Black spruce", "Other species")))[
                                                      ,':='(Year = -Inf, Yearend = -Inf, 
                                                            Value = Inf, Valueend = -Inf)]
treewaylinetype1 <- fixedEffect[Variable == "d1",][,':='(Direction = "IntraH", lineType = 1)]
treewaylinetype1[Lower95*Upper<=0, lineType:=2]
treewaylinetype2 <- fixedEffect[Variable == "d2",][,':='(Direction = "InterH", lineType = 1)]
treewaylinetype2[Lower95*Upper<=0, lineType:=2]
threewaylinetype <- rbind(treewaylinetype1, treewaylinetype2)[,.(Species, Direction, lineType)]
Figure3Data <- setkey(Figure3Data, Species, Direction)[setkey(threewaylinetype, Species, Direction),
                                                       nomatch = 0]
rm(treewaylinetype1, treewaylinetype2)
Figure3Data[, ':='(Direction = factor(Direction, levels = c("IntraH", "InterH")),
                   Species = factor(Species, levels = c("All species", "Jack pine", "Trembling aspen", 
                                                        "Black spruce", "Other species")),
                   lineType = factor(lineType, levels = c(1, 2)))]

Figure3 <- ggplot(data = Figure3Data[Main == 0,], aes(x = Year, y = Value))+
  geom_line(aes(col = RBI, group = RBI, linetype = lineType), size = 1)+
  geom_point(data = Figure3Data[Main == 1,], aes(x = Year, y = Value), col = "blue", size = 2)+
  geom_line(data = CompetitionTrends, aes(x = Year, y = Value, linetype = lineType), col = "blue", size = 1)+
  geom_errorbar(data = Figure3Data[Main == 1,], aes(ymin = Value_Lower, 
                                                    ymax = Value_Upper), 
                width = 1, col = "blue", size = 1)+
  guides(linetype = FALSE)+
  scale_colour_continuous(name = "RBI", low = "#FF0000", high = "#00FF00", breaks = c(0, 50, 100))+
  facet_grid(Species~Direction, scales = "free_y")+
  geom_segment(data = yaxisLineData,
               aes(x = Year, xend = Yearend, y = Value, yend = Valueend), size = 1.5)+
  geom_segment(aes(x = 1985, xend = 2010, y = 0, yend = 0), 
               col = "gray", linetype = 2, size = 1)+
  geom_text(data = data.frame(Direction = c("IntraH", "InterH"), labels = c("a", "b"), 
                              Species = "All species",
                              Value = -0.05, Year=1986),
            aes(x = Year, y = Value, label = labels), size = 10)+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  scale_y_continuous(name = "Effect of competition on tree growth", 
                     breaks = equal_breaks(n=4, s=0.05))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.background = element_rect(colour = "white"),
        strip.text.y = element_text(size = 16),
        strip.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.82, 0.85),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

dev(4)
clearPlot()

ggsave(file = file.path(workPath, "TablesFigures", "Figure 3_Bayesian_temporal trends in competition.png"),
       Figure3,  width = 10, height = 10)
