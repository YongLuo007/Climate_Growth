rm(list = ls())
library(SpaDES); library(ggplot2)
workPath <- "~/GitHub/Climate_Growth"
analysesData <- read.csv(file.path(workPath, "data", "MBdatafinal.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>% data.table
studyspecies <- c("JP", "TA", "BS")

# normality check
analysesData <- analysesData[Species %in% studyspecies, ]
# visual check
nomalityCheck <- ggplot(analysesData, aes(x = BiomassGR))+
  geom_histogram(aes(log(BiomassGR)), bins = 50)+
  facet_wrap(~Species)


# check the multilinearity among independent variables using Zuur 2010 method
source(file.path(workPath, "Rcodes", "Rfunctions", "HighstatLib_Zuur.R"))
for(indispecies in studyspecies){
  independentV1 <- data.frame(analysesData[Species == indispecies,
                                           .(RBI, Hegyi, IniDBH,
                                             Year)])

  viftest <- data.frame(corvif(independentV1))
  viftest$Species <- indispecies
  if(indispecies == "JP"){
    viftestoutput <- viftest
  } else {
    viftestoutput <- rbind(viftestoutput, viftest)
  }
}
print(viftestoutput)
GVIF Species
# RBI     1.386561      JP
# Hegyi   1.737999      JP
# IniDBH  1.653822      JP
# Year    1.053010      JP
# RBI1    1.350198      TA
# Hegyi1  2.049455      TA
# IniDBH1 2.216450      TA
# Year1   1.042794      TA
# RBI2    1.858667      BS
# Hegyi2  1.614244      BS
# IniDBH2 1.893845      BS
# Year2   1.028955      BS

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



