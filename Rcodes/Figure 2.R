rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)


workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "theBestModels.RData"))
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
analysesData <- analysesData[Species %in% c("JP", "TA", "BS"),]

ThreeDGrowthvsYearandDom <- data.table(Species = character(), Dominance = numeric(),
                                       Year = numeric(), PredictedBGR = numeric(),
                                       PredictedBGR_Upper = numeric(), PredictedBGR_Lower = numeric(),
                                       Main = numeric())
OtherVariable <- data.table(Species = character(), DBH = numeric(), H = numeric())
for(indispecies in c("JP", "BS", "TA")){
  themodel <- thebestmodel[[paste(indispecies, "_bestModel", sep = "")]]
  newdata <- data.table(expand.grid(Year = seq(min(analysesData[Species == indispecies,]$Year), 
                                               max(analysesData[Species == indispecies,]$Year),
                                               length = 100),
                                    Dominance = c(seq(0, 100, by = 5),
                                                  mean(analysesData[Species == indispecies,]$Dominance_indiBiomass))))
  newdata[,':='(Yearctd = Year-mean(analysesData[Species == indispecies,]$Year),
                Dominancectd = log(Dominance+1) - 
                  mean(log(analysesData[Species == indispecies,]$Dominance_indiBiomass+1)),
                logDBHctd = log(10)-
                  mean(log(analysesData[Species == indispecies,]$IniDBH)),
                logHctd = log(115)-
                  mean(log(analysesData[Species == indispecies,]$Hegyi)),
                Species = indispecies)]
  reducedFomu <- as.formula(paste(as.character(formula(summary(themodel)$term))[c(2, 1, 3)], collapse = " "))
  fittedvalues <- predict(themodel, newdata = newdata, level = 0, se.fit = TRUE)
  
  newdata$PredictedBGR <- exp(fittedvalues$fit)
  newdata$PredictedBGR_Upper <- exp(fittedvalues$fit+1.98*fittedvalues$se.fit)
  newdata$PredictedBGR_Lower <- exp(fittedvalues$fit-1.98*fittedvalues$se.fit)
  newdata[,Main:=0]
  newdata[Dominance == mean(analysesData[Species == indispecies,]$Dominance_indiBiomass), Main := 1]
  ThreeDGrowthvsYearandDom <- rbind(ThreeDGrowthvsYearandDom,
                                    newdata[,.(Species, Year, Dominance, PredictedBGR, 
                                               PredictedBGR_Upper, PredictedBGR_Lower, Main)])
  rm(themodel, newdata, fittedvalues, reducedFomu)
}


ThreeDGrowthvsYearandDom$Species <- factor(ThreeDGrowthvsYearandDom$Species,
                                           levels = c("JP", "TA", "BS"),
                                           labels = c("Jack pine", "Trembling aspen", "Black spruce"))
OtherVariable$Species <- factor(OtherVariable$Species,
                                levels = c("JP", "TA", "BS"),
                                labels = c("Jack pine", "Trembling aspen", "Black spruce"))
OtherVariable[, ':='(y1 = 5, y2 = 4, Year = 1993, DBH1 = paste(DBH, "cm"))]

Fig_a_3D <- ggplot(data = ThreeDGrowthvsYearandDom[Main == 0], aes(x = Year, y = PredictedBGR))+
  geom_line(aes(group = Dominance, col = Dominance))+
  scale_colour_continuous(low = "#FF0000", high = "#00FF00", breaks = seq(0, 100, by = 20))+
  geom_text(data = data.frame(Year = 1985, y = 3, label = "a", Species = "Jack pine"),
            aes(x = Year, y = y, label = label), size = 8)+
  geom_line(data = ThreeDGrowthvsYearandDom[Main == 1], aes(x = Year, y = PredictedBGR),
            colour = "black", size = 1)+
  facet_wrap(~Species)+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  scale_y_continuous(name = expression(atop("Aboveground biomass growth rate", paste("(Kg ", year^{-1}, ")"))),
                     limits = c(0, 3), breaks = seq(0, 3, by = 0.5))+
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
figure_domidensity <- ggplot(data = analysesData, aes(x = Dominance_indiBiomass))+
  geom_density(aes(col = Dominance_indiBiomass),fill = "red")+
  facet_wrap(~Species)



ggsave(file = file.path(workPath, "TablesFigures","simulatedBGRTemporalTrends.png"), Fig_a_3D,
       width = 13.3, height = 5.8, units = "in")



