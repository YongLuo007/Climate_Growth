rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "Results", "BiomassGR_mainResults.RData"))

# how to calculate the relative importance for each predictor

# figure 2 year effect and change with dominance
YearWithDominanceTable <- data.table(Species = character(), Dominance = numeric(),
                                    Value = numeric(), SE = numeric(), MainEffect = numeric())
ThreeDGrowthvsYearandDom <- data.table(Species = character(), Dominance = numeric(),
                                       Year = numeric(), PredictedBGR = numeric(),
                                       PredictedBGR_Upper = numeric(), PredictedBGR_Lower = numeric(),
                                       Main = numeric())
OtherVariable <- data.table(Species = character(), DBH = numeric(), H = numeric())
for(indispecies in c("JP", "BS", "TA")){
  set.seed(1)
  tempdata <- data.table(Dominance = seq(0, 100, by = 20))[
    , Dominctd:=Dominance-mean(analysesData[Species == indispecies,]$Dominance_Biomass)]
  yeareffect <- rnorm(10000,
                      mean = OverallResults[Species == indispecies & Variable == "Yearctd",]$Value,
                      sd = 100*OverallResults[Species == indispecies & Variable == "Yearctd",]$Std.Error)
  yearDominanceeffect <- rnorm(10000, 
                               mean = OverallResults[Species == indispecies & Variable == "Yearctd:Dominancectd",]$Value,
                               sd = 100*OverallResults[Species == indispecies & Variable == "Yearctd:Dominancectd",]$Std.Error)
  for(i in 1:nrow(tempdata)){
    tempdata$Value[i] <- mean(tempdata$Dominctd[i]*yearDominanceeffect+yeareffect)
    tempdata$SE[i] <- sd(tempdata$Dominctd[i]*yearDominanceeffect+yeareffect)/100
  }
  YearWithDominanceTable <- rbind(YearWithDominanceTable,
                                 tempdata[,.(Species=indispecies, Dominance, Value, SE, MainEffect = 0)],
                                 data.table(Species = indispecies, 
                                            Dominance = mean(analysesData[Species == indispecies,]$Dominance_Biomass),
                                            Value = OverallResults[Species == indispecies & Variable == "Yearctd",]$Value,
                                            SE = OverallResults[Species == indispecies & Variable == "Yearctd",]$Std.Error,
                                            MainEffect = 1))
  themodel <- OverallModels[[paste(indispecies, "_all", sep = "")]]
  newdata <- data.table(expand.grid(Year = seq(min(analysesData[Species == indispecies,]$Year), 
                                               max(analysesData[Species == indispecies,]$Year),
                                               length = 100),
                                    Dominance = c(seq(0, 100, by = 1),
                                                  mean(analysesData[Species == indispecies,]$Dominance_Biomass))))
  newdata[,':='(Yearctd = Year-mean(analysesData[Species == indispecies,]$Year),
                Dominancectd = Dominance - mean(analysesData[Species == indispecies,]$Dominance_Biomass),
                logDBHctd = log(mean(analysesData[Species == indispecies,]$IniDBH))-
                  mean(log(analysesData[Species == indispecies,]$IniDBH)),
                logHctd = log(mean(analysesData[Species == indispecies,]$Hegyi))-
                  mean(log(analysesData[Species == indispecies,]$Hegyi)),
                Species = indispecies)]
  fittedvalues <- predict(themodel, newdata = newdata, level = 0, se.fit = TRUE)
  
  newdata$PredictedBGR <- exp(fittedvalues$fit)
  newdata$PredictedBGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  newdata$PredictedBGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  newdata[,Main:=0]
  newdata[Dominance == mean(analysesData[Species == indispecies,]$Dominance_Biomass), Main := 1]
  ThreeDGrowthvsYearandDom <- rbind(ThreeDGrowthvsYearandDom,
                                    newdata[,.(Species, Year, Dominance, PredictedBGR, 
                                               PredictedBGR_Upper, PredictedBGR_Lower, Main)])
  OtherVariable <- rbind(OtherVariable, 
                         data.table(Species = indispecies,
                                    DBH = round(mean(analysesData[Species == indispecies,]$IniDBH), 2),
                                    H = round(mean(analysesData[Species == indispecies,]$Hegyi), 2)))
  rm(themodel, newdata, tempdata, yearDominanceeffect, yeareffect, indispecies)
}

YearWithDominanceTable[, ':='(Value = Value*10, SE = SE*10)]
pvalues <- OverallResults[Variable == "Yearctd:Dominancectd",.(Species, p.value)]
YearWithDominanceTable <- setkey(YearWithDominanceTable, Species)[setkey(pvalues, Species), nomatch = 0]
YearWithDominanceTable[,linetype := 1]
YearWithDominanceTable[p.value >= 0.05,linetype := 2]
YearWithDominanceTable$Species <- factor(YearWithDominanceTable$Species,
                                         levels = c("JP", "TA", "BS"),
                                         labels = c("Jack pine", "Trembling aspen", "Black spruce"))
YearWithDominanceTable$linetype <- factor(YearWithDominanceTable$linetype,
                                         levels = c(1, 2),
                                         labels = c(1, 2))

# the a
Fig_a <- ggplot(data = YearWithDominanceTable[MainEffect == 0,], 
                aes(x = Dominance, y = Value))+
  geom_segment(aes(x = 0, xend = 100, y = 0, yend = 0), 
               linetype = 2, size = 1, colour = "gray")+
  geom_ribbon(aes(ymin = Value-1.98*SE, ymax = Value+1.98*SE, 
                  fill = Species, group = Species), alpha = 0.1,
              show.legend = FALSE)+
  geom_line(aes(col = Species, group = Species), size = 1)+
  geom_point(data = YearWithDominanceTable[MainEffect == 1,],
             aes(x = Dominance, y = Value, col = Species), size = 2)+
  geom_errorbar(data = YearWithDominanceTable[MainEffect == 1,],
           aes(x = Dominance, ymin = Value-1.98*SE, ymax = Value+1.98*SE, col = Species))+
  scale_color_manual(name = "Species", 
                     values = c("red", "green", "blue"))+
  scale_y_continuous(name = expression(paste("Year effect on growth (", 10^{-2}, ")")),
                     limits = c(-1, 3), 
                     breaks = round(seq(-1, 3, by = 2)))+
  scale_x_continuous(name = "Dominance index", limits = c(0, 100), breaks = seq(0, 100, by = 20))+
  annotate("text", x = 1, y = 5, label = "a", size = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))
ggsave(file = file.path(workPath,  "TablesFigures", "YearEffect.png"), Fig_a, width = 10, height = 5, units = "in")

ThreeDGrowthvsYearandDom$Species <- factor(ThreeDGrowthvsYearandDom$Species,
                                           levels = c("JP", "TA", "BS"),
                                           labels = c("Jack pine", "Trembling aspen", "Black spruce"))
# ThreeDGrowthvsYearandDom[Main==1, Dominance:=0]

# ThreeDGrowthvsYearandDom$Dominance <- factor(ThreeDGrowthvsYearandDom$Dominance, levels = c(0, 0, 50, 100),
#                                              labels = c("DI = mean", "DI = 0", "DI = 50", "DI = 100"))

OtherVariable$Species <- factor(OtherVariable$Species,
                                levels = c("JP", "TA", "BS"),
                                labels = c("Jack pine", "Trembling aspen", "Black spruce"))

OtherVariable[, ':='(y1 = 5, y2 = 4, Year = 1993, DBH1 = paste(DBH, "cm"))]
Fig_a_3D <- ggplot(data = ThreeDGrowthvsYearandDom[Main == 0], aes(x = Year, y = PredictedBGR))+
  # geom_ribbon(aes(group = Dominance, fill = as.factor(Dominance), 
  #                 ymin = PredictedBGR_Lower, ymax = PredictedBGR_Upper), alpha = 0.1)+
  geom_line(aes(group = Dominance, col = Dominance))+
  scale_colour_continuous(low = "#FF0000", high = "#00FF00", breaks = seq(0, 100, by = 20))+
  # scale_colour_manual(name = "Dominance", values = rev(colorRamps::green2red(n=10)),
  #                     labels = seq(1, 10, by = 1))+
  # scale_fill_manual(name = "Dominance", values = rev(colorRamps::green2red(n=10)),
  #                   labels = seq(1, 10, by = 1))+
  geom_text(data = data.frame(Year = rep(1988, 6), y = c(rep(5, 3), rep(4, 3)), 
                              labels = c(rep("DBH:", 3), rep("H:", 3)),
                              Species = rep(c("Jack pine", "Trembling aspen", "Black spruce"), 2)),
            aes(x = Year, y = y, label = labels), hjust = 0, size = 5)+
  geom_text(data = OtherVariable, aes(x = Year, y = y1, label = DBH1), hjust = 0, size = 5)+
  geom_text(data = OtherVariable, aes(x = Year, y = y2, label = H), hjust = 0, size = 5)+
  geom_text(data = data.frame(Year = 1985, y = 8, label = "a", Species = "Jack pine"),
            aes(x = Year, y = y, label = label), size = 8)+
  facet_wrap(~Species)+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  scale_y_continuous(name = expression(atop("Aboveground biomass growth rate", paste("(Kg ", year^{-1}, ")"))),
                     limits = c(0, 10), breaks = seq(0, 6, by = 1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.background = element_rect(colour = "white"),
        strip.text = element_text(size = 15),
        legend.position = c(0.2, 0.75),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.direction = "horizontal")

ggsave(file = file.path(workPath, "TablesFigures","simulatedBGRTemporalTrends.png"), Fig_a_3D,
width = 13.3, height = 5.8, units = "in")


# figure 2 competition effect change with dominance
CompetitionWithDominanceTable <- data.table(Species = character(), 
                                           Dominance = numeric(),
                                    Value = numeric(), SE = numeric(), MainEffect = numeric())
for(indispecies in c("JP", "BS", "TA")){
  set.seed(1)
  tempdata <- data.table(Dominance = seq(0, 100, by = 10))[
    , Dominctd:=Dominance-mean(analysesData[Species == indispecies,]$DominanceIndex)]
  IntraHeffect <- rnorm(10000,
                        mean = OverallResults[Species == indispecies & Variable == "logHctd",]$Value,
                        sd = 100*OverallResults[Species == indispecies & Variable == "logHctd",]$Std.Error)
  IntraHDominanceeffect <- rnorm(10000, 
                                 mean = OverallResults[Species == indispecies & Variable == "logHctd:Dominancectd",]$Value,
                                 sd = 100*OverallResults[Species == indispecies & Variable == "logHctd:Dominancectd",]$Std.Error)
  for(i in 1:nrow(tempdata)){
    tempdata$Value[i] <- mean(tempdata$Dominctd[i]*IntraHDominanceeffect+IntraHeffect)
    tempdata$SE[i] <- sd(tempdata$Dominctd[i]*IntraHDominanceeffect+IntraHeffect)/100
  }
  CompetitionWithDominanceTable <- rbind(CompetitionWithDominanceTable,
                                 tempdata[,.(Species=indispecies, 
                                             Dominance, Value, SE, MainEffect = 0)],
                                 data.table(Species = indispecies,
                                            Dominance = mean(analysesData[Species == indispecies,]$DominanceIndex),
                                            Value = OverallResults[Species == indispecies & Variable == "logHctd",]$Value,
                                            SE = OverallResults[Species == indispecies & Variable == "logHctd",]$Std.Error,
                                            MainEffect = 1))
  
  
  
}
print( OverallResults[Variable == "logHctd:Dominancectd",.(Species, Value, p.value)])


Fig_b <- ggplot(data = CompetitionWithDominanceTable[MainEffect == 0,], 
                aes(x = Dominance, y = Value))+
  geom_segment(aes(x = 0, xend = 100, y = 0, yend = 0), colour = "gray",
               linetype = 2, size = 1)+
  geom_ribbon(aes(ymin = Value-1.98*SE, ymax = Value+1.98*SE, 
                  fill = as.factor(Species), group = Species), alpha = 0.1)+
  geom_line(aes(col = as.factor(Species), linetype = as.factor(Species), group = Species),
            size = 1)+
  geom_point(data = CompetitionWithDominanceTable[MainEffect == 1,],
             aes(x = Dominance, y = Value, col = as.factor(Species)), size = 2)+
  geom_errorbar(data = CompetitionWithDominanceTable[MainEffect == 1,],
                aes(x = Dominance, ymin = Value-1.98*SE, ymax = Value+1.98*SE, col = as.factor(Species)))+
  scale_color_manual(name = "Species", 
                     values = c("red", "green", "blue"),
                     labels = c("Black spruce", "Jack pine", "Trembling aspen"))+
  # OverallResults[Variable == "logIntraHctd",.(Species,Value, p.value)]
  scale_linetype_manual(name = "Species", 
                     values = c(1, 1, 1),
                     labels = c("Black spruce", "Jack pine", "Trembling aspen"))+
  scale_fill_manual(name = "Species", 
                    values = c("red", "green", "blue"),
                    labels = c("Black spruce", "Jack pine", "Trembling aspen"))+
  scale_y_continuous(name = "Effect of competition on growth",
                     limits = c(-1.5, 0.6), 
                     breaks = round(seq(-1.5, 0.6, by = 0.3), 1))+
  scale_x_continuous(name = "Dominance index", limits = c(0, 100), breaks = seq(0, 100, by = 20))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))




ggsave(file = file.path(workPath, "CompetitionEffectsWithDominance.png"), Fig_b, width = 10, height = 5)




# competition effect change with year

# figure 2 competition effect change with dominance
CompetitionWithYearTable <- data.table(Species = character(), 
                                            Year = numeric(),
                                            Value = numeric(), SE = numeric(), MainEffect = numeric())
for(indispecies in c("JP", "BS", "TA")){
  set.seed(1)
  YearRange <- range(analysesData[Species == indispecies,]$Year)
  tempdata <- data.table(Year = seq(YearRange[1], YearRange[2], length = 10))[
    , Yearctd:=Year-mean(analysesData[Species == indispecies,]$Year)]
  IntraHeffect <- rnorm(10000,
                        mean = OverallResults[Species == indispecies & Variable == "logHctd",]$Value,
                        sd = 100*OverallResults[Species == indispecies & Variable == "logHctd",]$Std.Error)
  IntraHDominanceeffect <- rnorm(10000, 
                                 mean = OverallResults[Species == indispecies & Variable == "Yearctd:logHctd",]$Value,
                                 sd = 100*OverallResults[Species == indispecies & Variable == "Yearctd:logHctd",]$Std.Error)
  for(i in 1:nrow(tempdata)){
    tempdata$Value[i] <- mean(tempdata$Yearctd[i]*IntraHDominanceeffect+IntraHeffect)
    tempdata$SE[i] <- sd(tempdata$Yearctd[i]*IntraHDominanceeffect+IntraHeffect)/100
  }
  CompetitionWithYearTable <- rbind(CompetitionWithYearTable,
                                         tempdata[,.(Species=indispecies, 
                                                     Year, Value, SE, MainEffect = 0)],
                                         data.table(Species = indispecies, 
                                                    Year = mean(analysesData[Species == indispecies,]$Year),
                                                    Value = OverallResults[Species == indispecies & Variable == "logHctd",]$Value,
                                                    SE = OverallResults[Species == indispecies & Variable == "logHctd",]$Std.Error,
                                                    MainEffect = 1))
  
  
}


OverallResults[Variable == "Yearctd:logHctd",.(Species, Value, p.value)]
# Species         Value      p.value
# 1:      JP -0.0134080947 3.874336e-15
# 2:      BS  0.0004601169 7.701253e-01
# 3:      TA  0.0073422123 2.362666e-03

Fig_c <- ggplot(data = CompetitionWithYearTable[MainEffect == 0 ,], 
                aes(x = Year, y = Value))+
  geom_ribbon(aes(ymin = Value-1.98*SE, ymax = Value+1.98*SE, 
                  fill = as.factor(Species), group = Species), alpha = 0.1)+
  geom_line(aes(col = as.factor(Species), linetype = as.factor(Species), group = Species),
            size = 1)+
  geom_point(data = CompetitionWithYearTable[MainEffect == 1,],
             aes(x = Year, y = Value, col = as.factor(Species)),
             size = 2)+
  geom_errorbar(data = CompetitionWithYearTable[MainEffect == 1,],
                aes(x = Year, ymin = Value-1.98*SE, ymax = Value+1.98*SE, col = as.factor(Species)))+
  scale_color_manual(name = "Species", 
                     values = c("red", "green", "blue"),
                     labels = c("Black spruce", "Jack pine", "Trembling aspen"))+
  # OverallResults[Variable == "YearbyIntraH",.(Species,Value, p.value)]
  scale_linetype_manual(name = "Species", 
                        values = c(2, 1, 1),
                        labels = c("Black spruce", "Jack pine", "Trembling aspen"))+
  scale_fill_manual(name = "Species", 
                    values = c("red", "green", "blue"),
                    labels = c("Black spruce", "Jack pine", "Trembling aspen"))+
  scale_y_continuous(name = "Effect of competition", 
                     limits = c(-0.94, -0.2), 
                     breaks = round(seq(-0.9, -0.2, by = 0.1), 1))+
  scale_x_continuous(name = "Year", limits = c(1987, 2008.5), breaks = seq(1987, 2010, by = 3))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = c(0.2, 0.2),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))



ggsave(file = file.path(workPath, "CompetitionEffectsWithYear.png"), Fig_c, width = 10, height = 5)






# the competition change with year by dominance
# three dimentional figure
CompetitionWithYearByDominanceTable <- data.table(Species = character(),
                                       Year = numeric(), Dominance = numeric(),
                                       Heffect = numeric())
Growth_H_ByYearDom <- data.table(Species = character(),
                                 Year = numeric(), Dominance = numeric(),
                                 Hegyi = numeric(), PredictedBGR = numeric())

for(indispecies in c("JP", "BS", "TA")){
  YearRange <- range(analysesData[Species == indispecies,]$Year)
  tempdata <- data.table(expand.grid(Year = seq(YearRange[1], YearRange[2], length = 100),
                                     Dominance = c(seq(0, 100, length = 100))))[
    , ':='(Yearctd = Year - mean(analysesData[Species == indispecies,]$Year),
           Domictd = Dominance - mean(analysesData[Species == indispecies,]$DominanceIndex))]
  MainHeffect <- OverallResults[Species == indispecies &
                              Variable == "logHctd",]$Value
  HChangewithYear <- OverallResults[Species == indispecies &
                                     Variable == "Yearctd:logHctd",]$Value
  HChangewithDominance <- OverallResults[Species == indispecies &
                                         Variable == "Yearctd:Dominancectd",]$Value
  HChangewithDomBYYear <- OverallResults[Species == indispecies &
                                                   Variable == "Yearctd:logHctd:Dominancectd",]$Value
  tempdata[,':='(Heffect = MainHeffect+Yearctd*HChangewithYear+Domictd*HChangewithDominance+
                                Yearctd*Domictd*HChangewithDomBYYear)]
  CompetitionWithYearByDominanceTable <- rbind(CompetitionWithYearByDominanceTable,
                                    tempdata[,.(Species=indispecies, 
                                                Year, Dominance, Heffect)])
  tempdata2 <- data.table(expand.grid(Year = c(1988, 1998, 2008),
                                     Dominance = c(0, 50, 100),
                                     Hegyi = seq(min(analysesData[Species == indispecies,]$Hegyi),
                                                 max(analysesData[Species == indispecies,]$Hegyi), 
                                                 length = 10000)))
  tempdata2[, ':='(Yearctd = Year - mean(analysesData[Species == indispecies,]$Year),
                   Dominancectd = Dominance - mean(analysesData[Species == indispecies,]$DominanceIndex),
                   logHctd = log(Hegyi)-mean(log(analysesData[Species == indispecies,]$Hegyi)),
                   logDBHctd = log(mean(analysesData[Species == indispecies,]$IniDBH))-
                     mean(log(analysesData[Species == indispecies,]$IniDBH)))]
  themodel <- OverallModels[[paste(indispecies, "_all", sep = "")]]
  tempdata2$PredictedBGR <- exp(predict(themodel, newdata = tempdata2, level = 0))
  Growth_H_ByYearDom <- rbind(Growth_H_ByYearDom, tempdata2[,.(Species = indispecies,
                                                               Year, Dominance, Hegyi, PredictedBGR)])
  rm(tempdata2, tempdata, HChangewithDomBYYear, HChangewithDominance, MainHeffect, HChangewithYear)
}

print(OverallResults[Variable == "Yearctd:logHctd:Dominancectd",.(Species,Value, p.value)])
# Species        Value      p.value
# 1:      JP 0.0008110226 4.563543e-27
# 2:      BS 0.0003944429 1.359574e-06
# 3:      TA 0.0003404419 2.517457e-04

CompetitionWithYearByDominanceTable$Species <- factor(CompetitionWithYearByDominanceTable$Species,
                                                      levels = c("JP", "TA", "BS"),
                                                      labels = c("Jack pine", "Trembling aspen", "Black spruce"))
Fig_d_3panel <- ggplot(data = CompetitionWithYearByDominanceTable,
                       aes(x = Year, y = Heffect))+
  geom_line(aes(group = Dominance, col = Dominance), size = 1)+
  scale_colour_continuous(low = "#FF0000", high = "#00FF00")+
  geom_text(data = data.frame(Year = 1985, y = -0.2, label = "b", Species = "Jack pine"),
            aes(x = Year, y = y, label = label), size = 8)+
  facet_wrap(~Species)+
  scale_y_continuous(name = "Effect of competition", 
                     limits = c(-1.001, -0.2), 
                     breaks = round(seq(-1, -0.2, by = 0.20), 1))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1), 
        axis.line.y = element_line(colour = "black", size = 1), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")


Growth_H_ByYearDom[,':='(Species = factor(Species, levels = c("JP", "TA", "BS"),
                                          labels = c("Jack pine", "Trembling aspen", "Black spruce")),
                         Year = factor(Year, levels = c(1988, 1998, 2008), 
                                       labels = c("Year = 1988", "Year = 1998", "Year = 2008")),
                         Dominance = factor(Dominance, levels = c(0, 50, 100),
                                            labels = c("Dominance = 0", "Dominance = 50",
                                                       "Dominance = 100")))]
options(scipen=10000)
Fig_d_9panel <- ggplot(data = Growth_H_ByYearDom, aes(x = Hegyi, y = PredictedBGR))+
  geom_line(aes(group = Year, col = Year))+
  facet_grid(Dominance ~ Species)+
  scale_x_log10(name = "Hegyi", breaks = c(1, 10, 100, 1000))+
  scale_y_log10(name = expression(atop("Aboveground biomass growth rate", paste("(Kg ", year^{-1}, ")"))),
                breaks = round(c(0.01, 1, 100), 2))+
  scale_colour_manual(name = "Year", values = c("red", "green", "blue"), 
                      labels = c("1988", "1998", 2008))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.background = element_rect(colour = "white"),
        strip.text = element_text(size = 15),
        legend.position = c(0.7, 0.07),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.direction = "horizontal")
  
Figd_Grob <- ggplotGrob(Fig_d)
Fige_Grob <- ggplotGrob(Fig_e)
Figf_Grob <- ggplotGrob(Fig_f)


Figd_Grob$widths <- Fige_Grob$widths
Figd_Grob$heights <- Figf_Grob$heights
Figf_Grob$widths <- Fige_Grob$widths
Fige_Grob$heights <- Figf_Grob$heights


dev(4)
clearPlot()
plotlayout <- rbind(c(1), c(2), c(3))
c <- grid.arrange(Figd_Grob, Fige_Grob, Figf_Grob,
                  layout_matrix = plotlayout)

ggsave(file = file.path(workPath, "CompetitionEffectsWithYearByDominance.png"), c,
       width = 10, height = 15)
setwd(orgPath)

