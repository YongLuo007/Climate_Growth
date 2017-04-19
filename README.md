# Competition mediates climate change-tree growth relationship
Examine the role of competition in expressing climate change-tree growth relationships

## Data

You will need to get the orginal data from Yong Luo. The data were not put into repo for two reasons: 1. data size is too large to upload; 2. restriction of data usage

## Data preparation
The data were selected using 

```r
source(file.path(".", "Rcodes", "MBdataselection.R"))
source(file.path(".", "Rcodes", "IDV preparation.R"))
source(file.path(".", "Rcodes", "assignCI.R"))
```


The Hegyi index were derived from  `HeghyiCICalculationModified.R` in Rfunctions folder.
The climate anomalies that derived from BioSIM is in data folder.

## Data analyses

### Examine temporal trends of tree growth
For each species group among all trees, Jack pine, trembling aspen, black spruce and minor species group, mixed models that account for random intercept at tree and plot level, and random slopes of Year at these levels had been developed using `lme` function in `nlme` package. The models were conducted as a full model forms.

```r
source(file.path(".", "Rcodes", "YearFullModels.R"))

```


### Examine climate changes and their associations with tree growth

Follow the same procedure as best Year models, the best climate models that used specific climate anomalies had been implemented.
```r
source(file.path(".", "Rcodes", "BestClimateModels.R"))
```
For detecting regrional climate changes, mixed models were used. The results and figure can be produced using following codes:
```r
source(file.path(".", "Rcodes", "Figure S4 regional climate changes.R"))
```

## Tables and figures presentations
To visualize study area and regional climate changes, the following codes were used.
```r
source(file.path(".", "Rcodes", "Figure 1 plot locations.R"))
```
To calculate the relative importance of predictors in best Year models,
```r
source(file.path(".", "Rcodes", "RelativeImportance.R"))
```
To produce temporal trends in tree growth
```r
source(file.path(".", "Rcodes", "Figure 3 temporal change in ABGR with competition.R"))
```
To present the climate change associations with tree growth,
```r
source(file.path(".", "Rcodes", "Figure 4 ClimateGrowthRelationships.R"))
```

