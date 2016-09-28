# Climate-tree growth relationship
Examine how climate change affected tree growth

## Data

You will need to get the orginal data from Yong Luo. The data were not put into repo for two reasons: 1. data size is too large to upload; 2. restriction of data usage

## Data preparation
The data were selected using 

```r
source(file.patch(".", "Rcodes", "MBdataselection.R"))
```

The independent variables (IDVs) were incorporated into the data using

```r
source(file.path(".", "Rcodes", "IDV preparation.R"))
```
For Hegyi index, a CItable(in data folder) had been prepared using function `HeghyiCICalculation.R` in Rfunctions folder.
The climate anomalies is in data folder

## Data inspection
The data was inspected for the normality of depedendt variable (visually check), linearity between dependent variable and independent variables (visually check), and variance inflaction factor (VIF).


The inspection was processed using `dataInspection.R`. The VIF was checked using method proposed by Zuur, A.F., et al. (2010)

## Data analyses

### Examine temporal trends of tree growth
For each species among Jack pine, Trembling aspen and Black spruce, mixed models that account for both random slope and random intercept had been developed using `lme` function in `nlme` package.


For each species, a best model was determined based on smallest AIC and was considered as the final model that was reported in the results.

### Examine climate changes and their associations with tree growth


## Tables and figures presentations

