#' To generate the bug file for jags
#'
#' @param coeffs character, specify the terms in the model
#' 
#' @param model character, whether is the year model (yearModel) or climate model (climateModel)
#' 
#' @param modelName character, specify the model name, default is model, which will generate a model.bug file
#' 
#' @param path character, choose where the bug model should be written
#' 
#' @param predictTrends logical, where the model should include predictions
#' 
#'
#' @return no object returned
#' 
#' @importFrom data.table data.table ':='
#' @importFrom dplyr left_join '%>%' 
#'
#' @note no note
#'
#' @seealso no
#'
#' @example 
#' \dontrun{
#' library(data.table)
#' path <- file.path(".")
#' inits <- read.csv(file.path(".", "Results", "YearcoeffInits.csv"),
#'                   header = TRUE, stringsAsFactors = FALSE) %>% data.table
#' indispecies <- "All"
#' inits <- inits[species == indispecies, ]
#' coeffs <- inits[coeffs != "a", ]$coeffs[c(1, 9, 10, 13)]
#' bugGenerator(coeffs = coeffs, model = "climateModel", modelName = "test",
#'              path = path, predictTrends = TRUE)
#' 
#' }
#' rm(list = ls())

#'
#' @export
#' @docType methods
#' @rdname bugGenerator
#'
#' @author Yong Luo
#'
setGeneric("bugGenerator",
           function(coeffs,
                    model,
                    modelName,
                    path,
                    predictTrends) {
             standardGeneric("bugGenerator")
           })

#' @export
#' @rdname bugGenerator
setMethod(
  "bugGenerator",
  signature = c(coeffs = "character",
                model = "character",
                modelName = "character",
                path = "character",
                predictTrends = "logical"),
  definition = function(coeffs,
                        model,
                        modelName,
                        path,
                        predictTrends){
    if(model == "yearModel"){
      baseTable <- data.table(coeff = c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "c1", "c2", "c3", 
                                         "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", 
                                         "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21"),
                              term = c("b1*logDBHctd[i]", "b2*Yearctd[i]", "b3*logIntraHctd[i]", "b4*logInterHctd[i]", 
                                        "b5*RBIctd[i]", "b6*logSActd[i]", "b7*logSBctd[i]", "c1*logDBHctd[i]*Yearctd[i]", 
                                        "c2*logDBHctd[i]*logIntraHctd[i]", "c3*logDBHctd[i]*logInterHctd[i]", 
                                        "c4*logDBHctd[i]*RBIctd[i]", "c5*logDBHctd[i]*logSActd[i]", "c6*logDBHctd[i]*logSBctd[i]", 
                                        "c7*Yearctd[i]*logIntraHctd[i]", "c8*Yearctd[i]*logInterHctd[i]", 
                                        "c9*Yearctd[i]*RBIctd[i]", "c10*Yearctd[i]*logSActd[i]", "c11*Yearctd[i]*logSBctd[i]", 
                                        "c12*logIntraHctd[i]*logInterHctd[i]", "c13*logIntraHctd[i]*RBIctd[i]", 
                                        "c14*logIntraHctd[i]*logSActd[i]", "c15*logIntraHctd[i]*logSBctd[i]", 
                                        "c16*logInterHctd[i]*RBIctd[i]", "c17*logInterHctd[i]*logSActd[i]", 
                                        "c18*logInterHctd[i]*logSBctd[i]", "c19*RBIctd[i]*logSActd[i]", 
                                        "c20*RBIctd[i]*logSBctd[i]", "c21*logSActd[i]*logSBctd[i]"))
      baseTable <- baseTable[coeff %in% coeffs, ]
      equations <- paste(" logY[i] ~ a +", paste(baseTable$term, collapse = " + ", sep = ""),
                         "+ u[Plot[i]] + q[Tree[i]]")
      priors <- paste("", baseTable$coeff, "~ dnorm(1,0.000001)\n")
      predictTrends <- TRUE
      if(predictTrends){
        trendsPrediction <- " for(m in 1:50){
   predictBAGR1[m] <- exp(a+b2*YearPredict[m])
   predictIntraHeffect[m] <- b3+c7*YearPredict[m]
   predictInterHeffect[m] <- b4+c8*YearPredict[m]
  }
  for(n in 1:500){
   predictY2[n] <- exp(a+b2*YearRBI_Yearctd[n]+b5*YearRBI_RBIctd[n]+c9*YearRBI_RBIctd[n]*YearRBI_Yearctd[n])
  }\n"
} else {
  trendsPrediction <- ""
    }
      cat("model
{ for (i in 1:NofObs){
  Y[i] ~ dlnorm(logY[i], tau_obs) \n",
  equations,"\n", " } \n",
  " for (j in 1:NofPlot){
   u[j] ~ dnorm(0, tau_plot_intercept)
  }
  for(k in 1:NofTree){
   q[k] ~ dnorm(0, tau_tree)
  }
  tau_obs ~ dgamma(0.001,0.001)
  tau_plot_intercept ~ dgamma(0.001,0.001)
  tau_tree ~ dgamma(0.001,0.001)
  a ~ dnorm(1,0.000001)\n",
  priors,
  trendsPrediction,
  "}",
  file = file.path(path, paste(modelName, ".bug", sep = "")))
      
} else if(model == "climateModel"){
  baseTable <- data.table(coeff = c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "c1", "c2", "c3", 
                                    "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", 
                                    "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21"),
                          term = c("b1*logDBHctd[i]", "b2*Climatectd[i]", "b3*logIntraHctd[i]", "b4*logInterHctd[i]", 
                                   "b5*RBIctd[i]", "b6*logSActd[i]", "b7*logSBctd[i]", "c1*logDBHctd[i]*Climatectd[i]", 
                                   "c2*logDBHctd[i]*logIntraHctd[i]", "c3*logDBHctd[i]*logInterHctd[i]", 
                                   "c4*logDBHctd[i]*RBIctd[i]", "c5*logDBHctd[i]*logSActd[i]", "c6*logDBHctd[i]*logSBctd[i]", 
                                   "c7*Climatectd[i]*logIntraHctd[i]", "c8*Climatectd[i]*logInterHctd[i]", 
                                   "c9*Climatectd[i]*RBIctd[i]", "c10*Climatectd[i]*logSActd[i]", "c11*Climatectd[i]*logSBctd[i]", 
                                   "c12*logIntraHctd[i]*logInterHctd[i]", "c13*logIntraHctd[i]*RBIctd[i]", 
                                   "c14*logIntraHctd[i]*logSActd[i]", "c15*logIntraHctd[i]*logSBctd[i]", 
                                   "c16*logInterHctd[i]*RBIctd[i]", "c17*logInterHctd[i]*logSActd[i]", 
                                   "c18*logInterHctd[i]*logSBctd[i]", "c19*RBIctd[i]*logSActd[i]", 
                                   "c20*RBIctd[i]*logSBctd[i]", "c21*logSActd[i]*logSBctd[i]"))
  baseTable <- baseTable[coeff %in% coeffs, ]
  equations <- paste(" logY[i] ~ a +", paste(baseTable$term, collapse = " + ", sep = ""),
                     "+ u[Plot[i]] + q[Tree[i]]")
  priors <- paste("", baseTable$coeff, "~ dnorm(1,0.000001)\n")
  predictTrends <- TRUE
  if(predictTrends){
    trendsPrediction <- " for(m in 1:50){
    predictBAGR1[m] <- exp(a+b2*ClimatePredict[m])
    predictIntraHeffect[m] <- b3+c7*ClimatePredict[m]
    predictInterHeffect[m] <- b4+c8*ClimatePredict[m]
  }
  for(n in 1:500){
    predictY2[n] <- exp(a+b2*ClimateRBI_Climatectd[n]+b5*ClimateRBI_RBIctd[n]+c9*ClimateRBI_RBIctd[n]*ClimateRBI_Climatectd[n])
    }\n"
} else {
  trendsPrediction <- ""
}
cat("model
{ for (i in 1:NofObs){
  Y[i] ~ dlnorm(logY[i], tau_obs) \n",
  equations,"\n", " } \n",
  " for (j in 1:NofPlot){
  u[j] ~ dnorm(0, tau_plot_intercept)
  }
  for(k in 1:NofTree){
  q[k] ~ dnorm(0, tau_tree)
  }
  tau_obs ~ dgamma(0.001,0.001)
  tau_plot_intercept ~ dgamma(0.001,0.001)
  tau_tree ~ dgamma(0.001,0.001)
  a ~ dnorm(1,0.000001)\n",
  priors,
  trendsPrediction,
  "}",
  file = file.path(path, paste(modelName, ".bug", sep = "")))
}
})


