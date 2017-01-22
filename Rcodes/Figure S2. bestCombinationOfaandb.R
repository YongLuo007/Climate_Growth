rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(MuMIn);library(parallel)

output <- read.csv(file.path("~/GitHub/Climate_Growth/Results/bestWeightsbothsizeanddistanceRescaled.csv"),
                   header = T, stringsAsFactors = F) %>% data.table
output <- output[Species != "Other species",]

a <- melt(output, id.vars = c("Species", "sizeWeight", "disweight"), 
          measure.vars = c("IntraHAIC", "InterHAIC"),
          value.name = "Value")
a[,':='(variable = factor(variable, levels = c("IntraHAIC", "InterHAIC"),
                          labels = c("Intraspecific competition", "Interspecific competition")),
        Species = factor(Species, levels = c("Jack pine", "Trembling aspen", "Black spruce")))]

a[,':='(minvalue=min(Value), maxvalue = max(Value)), by = c("Species", "variable")]
a[, scaledValue:=(Value-minvalue)/(maxvalue-minvalue)]
cutpoints <- c(0, seq(0.05, 0.95, by = 0.1), 1)

minvaluepoints <- a[Value == minvalue,]
rectTable <- data.table(expand.grid(Species = c("Jack pine", "Trembling aspen", "Black spruce"),
                                    variable = c("Intraspecific competition", "Interspecific competition")))
rectTable[,':='(disweight = 0, xend = 0.7, sizeWeight = 6.5, yend = 8)]
alphaTable <- data.table::copy(minvaluepoints)[,':='(disweight = 0.1, alpha = paste("alpha==", sizeWeight))]
alphaTable[, sizeWeight := 7.25]
betaTable <- data.table::copy(minvaluepoints)[,':='(disweight = 0.4, alpha = paste("beta==", disweight))]
betaTable[, sizeWeight := 7.25]
figure <- ggplot(data=a, aes(x = disweight, y = sizeWeight))+
  facet_grid(Species~variable, scales = "free_y")+
  geom_raster(aes(fill = scaledValue), interpolate = TRUE)+
  geom_contour(aes(z = scaledValue), col = "gray50")+
  geom_point(data = minvaluepoints, aes(x = disweight, 
                                        y = sizeWeight),
             show.legend = FALSE, size = 3, col = "red")+
  geom_point(aes(x = 1, y = 0), col = "green",size = 3)+
  geom_rect(data = rectTable, aes(xmin = disweight, xmax = xend, ymin = sizeWeight, ymax = yend), 
            fill = "white", col = "black")+
  geom_text(data = alphaTable, 
             aes(x = disweight, y = sizeWeight, label = alpha),
             hjust = 0, size = 6, parse = TRUE)+
  geom_text(data = betaTable, 
          aes(x = disweight, y = sizeWeight, label = alpha),
          hjust = 0, size = 6, parse = TRUE)+
  scale_y_continuous(name = expression(paste("Cross-forests assymetric competition coefficient (", alpha, ")")))+
  scale_x_continuous(name = expression(paste("Crowdedness competition coefficient (", beta, ")")))+
  scale_fill_continuous(name = "Scaled AIC", breaks = c(0, 1))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 20),
        strip.background = element_rect(colour = "white"))

 ggsave(file = file.path("~/GitHub/Climate_Growth/TablesFigures/bestdistanceweight.png"),
        figure, width = 14, height = 16)

