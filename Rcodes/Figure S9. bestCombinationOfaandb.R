rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(MuMIn);library(parallel)
workPath <- "~/GitHub/Climate_Growth"
selectionMethod <- "AllCensus_PositiveGrowth_RandomPlotADTree"
output <- fread(file.path(workPath, "data", selectionMethod,
                             "bestAandB.csv"))
output[,':='(sizeWeight = as.numeric(sizeWeight), disWeight = as.numeric(disWeight))]
output <- output[sizeWeight<=8,]
a <- melt(output, id.vars = c("Species", "sizeWeight", "disWeight"), 
          measure.vars = c("HAIC"),
          value.name = "Value")
a <- rbind(a[1,][,Species := ""], a)
a[Species == "All species", Species:="All trees"]
a[,':='(Species = factor(Species, levels = c("All trees", "", "Jack pine",
                                             "Trembling aspen", "Black spruce",
                                             "Minor species")))]

a[,':='(minvalue=min(Value), maxvalue = max(Value)), by = c("Species", "variable")]
a[, scaledValue:=(Value-minvalue)/(maxvalue-minvalue)]
a[Species == "", scaledValue:="0"]
a[,sizeWeight:=as.numeric(sizeWeight)]
minvaluepoints <- a[Value == minvalue,][,':='(sizeWeight=as.numeric(sizeWeight))]
rectTable <- data.table(Species = factor(c("All trees", "Jack pine", 
                                                "Trembling aspen", "Black spruce",
                                                "Minor species"),
                                         levels = c("All trees", "", "Jack pine",
                                                    "Trembling aspen", "Black spruce",
                                                    "Minor species")))
rectTable[,':='(disWeight = 0, xend = 0.9, sizeWeight = 6.5, yend = 8)]
alphaTable <- data.table::copy(minvaluepoints)[,':='(disWeight = 0.1, alpha = paste("alpha==", sizeWeight))]
alphaTable[, sizeWeight := 7.25]
betaTable <- data.table::copy(minvaluepoints)[,':='(disWeight = 0.5, alpha = paste("beta==", disWeight))]
betaTable[, sizeWeight := 7.25]


figure <- ggplot(data=a[Species != "",], aes(x = disWeight, y = sizeWeight))+
  facet_wrap(~Species, scales = "free_y", ncol = 2, drop = FALSE)+
  geom_raster(aes(fill = scaledValue), interpolate = TRUE)+
  geom_contour(aes(z = scaledValue), col = "gray50")+
  geom_point(data = minvaluepoints[Species != "",], aes(x = disWeight, 
                                        y = sizeWeight),
             show.legend = FALSE, size = 3, col = "red")+
  geom_point(aes(x = 1, y = 0), col = "green",size = 3)+
  geom_rect(data = rectTable, aes(xmin = disWeight, xmax = xend, ymin = sizeWeight, ymax = yend), 
            fill = "white", col = "black")+
  geom_text(data = alphaTable[Species != "",], 
             aes(x = disWeight, y = sizeWeight, label = alpha),
             hjust = 0, size = 6, parse = TRUE)+
  geom_text(data = betaTable[Species != "",], 
          aes(x = disWeight, y = sizeWeight, label = alpha),
          hjust = 0, size = 6, parse = TRUE)+
  scale_y_continuous(name = expression(paste("Cross-stands assymetric coefficient (", alpha, ")")))+
  scale_x_continuous(name = expression(paste("Crowdedness coefficient (", beta, ")")))+
  scale_fill_continuous(name = "AIC\n", breaks = c(0, 1), labels = c("Minimum\n", "Maximum"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.background = element_rect(colour = "black"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 20, hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"))

 ggsave(file = file.path("~/GitHub/Climate_Growth/TablesFigures/Figure S9. bestdistanceweight.png"),
        figure, width = 12, height = 14)

