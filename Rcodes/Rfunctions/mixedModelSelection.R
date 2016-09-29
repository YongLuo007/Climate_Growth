#' the function to select the best model based on AIC, DIC or BIC
#'
#'
#' @param DV character, specify the dependent variable
#' 
#' @param IDV character vector, specify the independent variables
#' 
#' @param maxInteraction numeric, specify the maximum interaction among predictors
#'                                default is 1, which means no interaction terms
#'
#' @param ICTerm character specify the which information critiria that will been involved in 
#'                         the selection. e.g., AIC, BIC, DIC
#'                         default is AIC. 
#'                         
#' 
#' @param ... other arguements in nlme::lme function, excluding formula.
#'                  e.g., data, random, control and so on
#'                            
#' 
#' 
#' 
#' 
#'
#' @return a data table that has three columns, i.e., active, mapcode and ecoregion
#' 
#' @importFrom data.table data.table ':='
#' @importFrom dplyr left_join '%>%' 
#' @importFrom nlme lme
#' @importFrom MuMIn DIC
#' @importFrom stats AIC BIC
#'
#' @note no note
#'
#' @seealso no
#'
#' @export
#' @docType methods
#' @rdname mixedModelSelection
#'
#' @author Yong Luo
#'
setGeneric("mixedModelSelection",
           function(DV,
                    IDV,
                    maxInteraction,
                    ICTerm,
                    ...) {
             standardGeneric("mixedModelSelection")
           })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                ICTerm = "character"),
          definition = function(DV,
                                IDV,
                                maxInteraction,
                                ICTerm,
                                ...){
            allIDV <- IDV
            output <- data.table(Model = character(),
                                 Formula = character(),
                                 IC = numeric(),
                                 deltaIC = numeric(),
                                 MarR2 = numeric(),
                                 ConR2 = numeric())
            if(maxInteraction > 1){
              for(i in 2:maxInteraction){
                tempV <- data.frame(t(combinat::combn(IDV, i)))
                interactions <- as.character(tempV[,"X1"])
                for(j in 2:i){
                  interactions <- c(paste(interactions, ":",
                                          as.character(tempV[,paste("X", j, sep = "")]),
                                          sep = ""))
                }
                allIDV <- c(allIDV, interactions)
              }
            }
            rm(tempV, i, j, interactions)
            # full model
            modelFormula <- paste(DV, "~", paste(allIDV, collapse = "+"), sep = "")
            themodel <- nlme::lme(as.formula(modelFormula),...)
            tTable <- data.frame(summary(themodel)$tTable)
            reducedIDV <- row.names(tTable)[tTable$p.value < 0.05 & row.names(tTable) != "(Intercept)"]
            outputAdd <- data.table(Model = "Full",
                                    Formula = "FULL",
                                    IC = getIC(model = themodel, x = ICTerm),
                                    deltaIC = 0,
                                    MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                                    ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
            output <- rbind(output, outputAdd)
            rm(tTable, outputAdd, modelFormula, themodel)
            prevvari <- as.character()
            for(i in 1:length(reducedIDV)){
              if(length(reducedIDV) != length(prevvari)){
                reducedFomu <- paste(DV, "~", paste(reducedIDV, collapse = "+"), sep = "")
                prevvari <- reducedIDV
                themodel <- nlme::lme(as.formula(reducedFomu),...)
                tTable <- data.frame(summary(themodel)$tTable)
                reducedIDV <- row.names(tTable)[tTable$p.value < 0.05 & row.names(tTable) != "(Intercept)"]
                if(length(reducedIDV) == length(prevvari)){
                  BaseIDV <- paste(reducedIDV, collapse = ", ")
                  BaseIC <- getIC(themodel, x = ICTerm)
                  BaseIDVindi <- reducedIDV
                  thebestmodel <- themodel
                  thebestIDV <- reducedIDV
                  outputAdd <- data.table(Model = "BaseModel",
                                          Formula = reducedFomu,
                                          IC = BaseIC,
                                          deltaIC = 0,
                                          MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                                          ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
                  output <- rbind(output, outputAdd)
                  output[, deltaIC:=IC-BaseIC]
                  rm(outputAdd)
                }
                rm(tTable, reducedFomu, themodel)
              }
            }
            
            rm(i, prevvari)
            # drop variables from base model and check the ICs
            preIDV <- character()
            for(i in 1:length(thebestIDV)){
              if(length(preIDV) != length(thebestIDV)){
                preIDV <- thebestIDV
                a <- dropOneVariableFun(output = output, testModel = thebestmodel,
                                      testIDV = thebestIDV, DV = DV, ICTerm = ICTerm, ...)
                thebestIDV <- a$thebestIDV
                thebestmodel <- a$thebestmodel
                output <- a$output
              }
            }
            rm(preIDV, a)
            preIDV <- character()
            for(i in 1:length(allIDV)){
              if(length(preIDV) != length(thebestIDV)){
                preIDV <- thebestIDV
                a <- addOneVariableFun(allIDV = allIDV, output = output, 
                                       testModel = thebestmodel,
                                        testIDV = thebestIDV, DV = DV, ICTerm = ICTerm, ...)
                thebestIDV <- a$thebestIDV
                thebestmodel <- a$thebestmodel
                output <- a$output
              }
            }
            output$deltaIC <- output$IC-output[Model == "BaseModel",]$IC
            setnames(output, c("IC", "deltaIC"), paste(c("", "delta"), ICTerm, sep = ""))
            return(invisible(list(modelSummary = output, bestModel = thebestmodel,
                                  bestIDV = thebestIDV)))
          })


getIC <- function(model, x){
  if(x == "AIC"){
    IC <- as.numeric(AIC(model))
  } else if (x == "BIC"){
    IC <- as.numeric(BIC(model))
  } else if (x == "DIC"){
    IC <- as.numeric(MuMIn::DIC(model))
  }
  return(IC)
}

dropOneVariableFun <- function(output, testModel, testIDV, DV, ICTerm, ...){
  BaseIC <- getIC(testModel, x = ICTerm)
  thebestIDV <- testIDV
  thebestmodel <- testModel
  if(length(testIDV) > 1){
    tempV <- data.frame(t(combinat::combn(testIDV,
                                          (length(testIDV)-1))), # drop one variable
                        stringsAsFactors = FALSE)
    for(i in 1:nrow(tempV)){
      reducedIDV <- as.character(tempV[i,])
      reducedFomu <- paste(DV, "~", paste(reducedIDV, collapse = "+"), sep = "")
      themodel <- nlme::lme(as.formula(reducedFomu),...)
      if(getIC(themodel, x = ICTerm)-BaseIC < 2 & getIC(themodel, x = ICTerm) < min(output$IC)){
        thebestmodel <- themodel
        thebestIDV <- reducedIDV
      }
      outputAdd <- data.table(Model = "ReducedModel",
                              Formula = reducedFomu,
                              IC = getIC(themodel, x = ICTerm),
                              deltaIC = getIC(themodel, x = ICTerm) - BaseIC,
                              MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                              ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
      output <- rbind(output, outputAdd)
      rm(outputAdd, reducedFomu, reducedIDV, themodel)
    }
    rm(i, tempV)
  } else {
    reducedFomu <- paste(DV, "~NULL", sep = "")
    themodel <- nlme::lme(as.formula(reducedFomu),...)
    outputAdd <- data.table(Model = "NULL",
                            Formula = reducedFomu,
                            IC = getIC(themodel, x = ICTerm),
                            deltaIC = getIC(themodel, x = ICTerm) - BaseIC,
                            MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                            ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
    output <- rbind(output, outputAdd)
    if(getIC(themodel, x = ICTerm)-BaseIC <= 2){
      thebestmodel <- themodel
      thebestIDV <- "NULL"
    }
    rm(reducedFomu, themodel, outputAdd)
  }
  return(list(thebestmodel = thebestmodel, thebestIDV = thebestIDV, output = output))
}


addOneVariableFun <- function(allIDV, output, testModel, testIDV, DV, ICTerm, ...){
  addedIDVs <- allIDV[!(allIDV %in% testIDV)]
  BaseIC <- getIC(testModel, x = ICTerm)
  thebestIDV <- testIDV
  thebestmodel <- testModel
  for(i in 1:length(addedIDVs)){
    addedIDV <- c(testIDV, addedIDVs[i])
    addedFomu <- paste(DV, "~", paste(addedIDV, collapse = "+"), sep = "")
    themodel <- nlme::lme(as.formula(addedFomu),...)
    if(getIC(themodel, x = ICTerm)-BaseIC < -2 & getIC(themodel, x = ICTerm) < min(output$IC)){
      thebestmodel <- themodel
      thebestIDV <- reducedIDV
    }
    outputAdd <- data.table(Model = "ExpandedModel",
                            Formula = addedFomu,
                            IC = getIC(themodel, x = ICTerm),
                            deltaIC = getIC(themodel, x = ICTerm) - BaseIC,
                            MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                            ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
    output <- rbind(output, outputAdd)
    rm(outputAdd, addedFomu, addedIDV, themodel)
  }
  rm(i)
  return(list(thebestmodel = thebestmodel, thebestIDV = thebestIDV, output = output))
}






#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                ICTerm = "character"),
          definition = function(DV,
                                IDV,
                                ICTerm,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction = 1, ICTerm, ...))
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                ICTerm = "missing"),
          definition = function(DV,
                                IDV,
                                maxInteraction,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction, ICTerm = "AIC", ...))
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                ICTerm = "missing"),
          definition = function(DV,
                                IDV,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction = 1, ICTerm = "AIC", ...))
          })