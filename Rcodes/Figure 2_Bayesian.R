rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
load(file.path(workPath, "data", "Bayesian_YearModels.RData"))


Figure2Data[,Species:=factor(Species, levels = c("All", "JP",
                                                 "TA", "BS",
                                                 "Other"),
                               labels = c("All species", "Jack pine",
                                          "Trembling aspen", "Black spruce",
                                          "Other species"))]
Figure2Data[, Direction:=factor(Direction, 
                                  levels = c("Overall Trend", "Change with RBI"),
                                  labels = c("Overall trend", "Change with RBI"))]
# Yeareffect[, Species := factor(Species, levels = c("All species", "Jack pine",
#                                                    "Trembling aspen", "Black spruce",
#                                                    "Other species"),
#                                labels = c("All species", "Jack pine",
#                                           "Trembling aspen", "Black spruce",
#                                           "Other species"))]
# YeareffectWithRBI[, Species := factor(Species, levels = c("All species", "Jack pine",
#                                                           "Trembling aspen", "Black spruce",
#                                                           "Other species"),
#                                       labels = c("All species", "Jack pine",
#                                                  "Trembling aspen", "Black spruce",
#                                                  "Other species"))]

AFCtable <- fixedEffect[Variable == "Yearctd",][,.(Species = factor(Species, levels = c("All", "JP",
                                                                                        "TA", "BS",
                                                                                        "Other"),
                                                                    labels = c("All species", "Jack pine",
                                                                                   "Trembling aspen", "Black spruce",
                                                                                   "Other species")), 
                                              Direction = factor("Overall trend", levels = c("Overall trend", "Change with RBI")),
                                              label = paste(round(exp(Mean)-1, 3),
                                                            " (", round(exp(Lower95)-1, 3), "-",
                                                            round(exp(Upper95)-1, 3), ")", sep = ""))]
AFCtable <- rbind(data.table::copy(AFCtable)[, ':='(Year = 1995, y = 0.95)], 
                  data.table::copy(AFCtable)[,':='(label = "AFC: ", Year = 1990, y = 0.95)])
Figure2 <- ggplot(data = Figure2Data, aes(x = Year, y = PredictedABGR))+
  # for a overall trends
  geom_ribbon(data = Figure2Data[Direction == "Overall trend",],
              aes(x = Year, ymin = PredictedABGR_Lower, ymax = PredictedABGR_Upper),
              fill = "blue", alpha = 0.1)+
  geom_line(data = Figure2Data[Direction == "Overall trend",],
            aes(x = Year, y = PredictedABGR), size = 1, col = "blue")+
  # for b change with RBI
  geom_line(data = Figure2Data[Direction == "Change with RBI",],
            aes(x = Year, y = PredictedABGR, group = RBI, col = RBI),
            size = 1)+
  # accessories
  geom_segment(aes(x=-Inf, xend=-Inf, y=Inf, yend=-Inf), size = 1.5)+
  geom_text(data = data.frame(Year = rep(1986, 2), y = 1.05, label = c("a", "b"), 
                              Species = "All species", Direction = c("Overall trend", "Change with RBI")),
            aes(x = Year, y = y, label = label), size = 10)+
  geom_text(data = AFCtable, aes(x = Year, y = y, label = label), hjust = 0, size = 5)+
  scale_colour_continuous(name = "RBI", low = "#FF0000", high = "#00FF00", breaks = c(0, 50, 100))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  scale_y_continuous(name = expression(paste("Aboveground biomass growth rate (Kg  ", year^{-1}, ")")),
                     limits = c(0, 1.07), breaks = seq(0.1, 1, 0.3))+
  
  facet_grid(Species~Direction)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        # axis.line.y = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.background = element_rect(colour = "white"),
        strip.text.y = element_text(size = 16),
        # strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.82, 0.78),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

dev(4)
clearPlot()
ggsave(file = file.path(workPath, "TablesFigures",
                        "Figure 2_Bayesian_temporal trends in ABGR.png"), 
       Figure2,
       width = 10, height = 10)
