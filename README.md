# Competition mediates climate change-tree growth relationship
Examine the role of competition in expressing climate change-tree growth relationship

## Data

You will need to get the orginal data from Yong Luo. The data were not put into repo for two reasons: 1. data size is too large to upload; 2. restriction of data usage

## Data preparation
The data were selected using 

```r
source(file.path(".", "Rcodes", "MBdataselection.R"))
```

The independent variables (IDVs) were incorporated into the data using

```r
source(file.path(".", "Rcodes", "IDV preparation.R"))
```
The Hegyi index were derived from  `HeghyiCICalculationModified.R` in Rfunctions folder.
The climate anomalies that derived from BioSIM is in data folder.

## Data analyses

### Examine temporal trends of tree growth
For each species among Jack pine, Trembling aspen and Black spruce, mixed models that account for random intercept had been developed using `lme` function in `nlme` package.
In the model selection process, the backward model selection based on AIC approach was used. The process was implemented using `mixedModelSelection.R`

```r
source(file.path(".", "Rcodes", "BestYearModels.R"))
```
For each species, a best model was determined based on smallest AIC and was considered as the final model that was reported in the results.

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

