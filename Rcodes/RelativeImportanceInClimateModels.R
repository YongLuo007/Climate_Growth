rm(list = ls())
# produce figure 2
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(gridExtra)
workPath <- "~/GitHub/Climate_Growth"
figuredata <- read.csv(file.path(workPath, "Results", "RelativeImportanceClimateModel.csv"),
                       header = TRUE, stringsAsFactors = FALSE) %>% data.table
figuredata[,':='(Species = unlist(lapply(Species_Climate, function(s) unlist(strsplit(s, "_", fixed = TRUE))[1])),
                 Climate = unlist(lapply(Species_Climate, function(s) unlist(strsplit(s, "_", fixed = TRUE))[2])))]
figuredata[,':='(Species = factor(Species,
                                  levels = c("Jack pine", 
                                             "Trembling aspen", "Black spruce", 
                                             "Other species")),
                 Group = factor(Group,levels = c("Ontogeny", "Competition",
                                                 "Ontogeny+Competition", "Climate", 
                                                 "Competition+Climate", "Ontogeny+Climate"),
                                labels = c("Ontogeny", "Competition",
                                           "Ontogeny×Competition", "Climate", 
                                           "Climate×Competition",
                                           "Climate×Ontogeny")))]
temperature <- c("ATA", "GSTA", "NONGSTA")
CMI <- c("ACMIA", "GSCMIA")
CO2 <- c("ACO2A")
figuredata1 <- figuredata[Climate == "ATA",]
fig1 <- ggplot(data = figuredata1, aes(x = Group, y = importance))+
  facet_grid(Species~., nrow = 2)+
  geom_segment(aes(xend = Group, yend = 0), size = 5)
