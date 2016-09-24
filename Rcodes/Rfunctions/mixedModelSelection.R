#' the function to select the best model based on AIC, DIC or BIC
#'
#' @param data data.table, the data
#'
#' @param DV character, specify the dependent variable
#' 
#' @param IDV character vector, specify the independent variables
#' 
#' @param maxInteraction numeric, specify the maximum interaction among predictors
#'                                default is 1, which means no interaction terms
#' 
#' @param randomMethod character, choose the random slope model by using "slope",
#'                            choose the random intercept model by using "intercept",
#'                            choose the both random slope and intercept model by using "both"
#'                            default is "intercept"
#'                            
#' @param randomTerm character, specify the higher hierachical layers
#' 
#' @param slopeTerm character, if randomMethod is set as "slope" or "both", a slopeTerm 
#'                             must be needed
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
           function(data,
                    DV,
                    IDV,
                    maxInteraction,
                    randomMethod,
                    randomTerm,
                    slopeTerm) {
             standardGeneric("mixedModelSelection")
           })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(data = "data.table",
                                DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                randomMethod = "character",
                                randomTerm = "character",
                                slopeTerm = "character"),
          definition = function(data,
                                DV,
                                IDV,
                                maxInteraction,
                                randomMethod,
                                randomTerm,
                                slopeTerm){
            allIDV <- IDV
            output <- data.table(Model = character(), Formula = character(),
                                 DIC = numeric(), AIC = numeric(), BIC = numeric(),
                                 MarR2 = numeric(), ConR2 = numeric())
            outputModel <- list()
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
            modelFormula <- paste(DV, "~", allIDV[1], sep = "")
            if(length(allIDV)>=2){
              for(j in 2:length(allIDV)){
                modelFormula <- paste(modelFormula, "+", allIDV[j], sep = "")
              }
            }
            if(randomMethod == "slope"){
              randomArg <- paste("random = ~", slopeTerm, "|", randomTerm, sep = "")
            } else if (randomMethod == "intercept"){
              randomArg <- paste("random = ~1|", randomTerm, sep = "")
            } else if (randomMethod == "both"){
              randomArg <- paste("random = ~1+", slopeTerm, "|", randomTerm, sep = "")
            }
            browser()
            themodelForm <- paste("themodel <- lme(", modelFormula, ",",
                                  randomArg, ",",
                                  "data = data,
                                  control = lmeControl(maxIter=10000, msMaxIter = 10000))", sep = "")
            eval(parse(text=themodelForm))
            outputModel[["Full"]] <- themodel
            rm(themodelForm)
            outputAdd <- data.table(Model = "Full",
                                    Formula = modelFormula,
                                    DIC = as.numeric(DIC(themodel)),
                                    AIC = as.numeric(AIC(themodel)),
                                    BIC = as.numeric(BIC(themodel)),
                                    MarR2 = as.numeric(r.squaredGLMM(themodel)[1]),
                                    ConR2 = as.numeric(r.squaredGLMM(themodel)[2]))
            output <- rbind(output, outputAdd)
            tTable <- data.frame(summary(themodel)$tTable)
            reducedIDV <- row.names(tTable)[tTable$p.value<0.05 & row.names(tTable) != "(Intercept)"]
            rm(tTable, outputAdd, j, modelFormula, themodel)
            prevvari <- as.character()
            for(i in 1:length(reducedIDV)){
              if(length(reducedIDV) != length(prevvari)){
                reducedFomu <- paste(DV, "~", reducedIDV[1], sep = "")
                if(length(reducedIDV)>=2){
                  for(j in 2:length(reducedIDV)){
                    reducedFomu <- paste(reducedFomu, "+", reducedIDV[j], sep = "")
                    rm(j)
                  }
                }
                themodelForm <- paste("themodel <- lme(", reducedFomu, ",",
                                      randomArg, ",",
                                      "data = data,
                                      control = lmeControl(maxIter=10000, msMaxIter = 10000))", sep = "")
                eval(parse(text=themodelForm))
                outputModel[[paste("ReducedModel", i, sep = "")]] <- themodel
                rm(themodelForm)
                outputAdd <- data.table(Model = paste("ReducedModel", i, sep = ""),
                                        Formula = reducedFomu,
                                        DIC = as.numeric(DIC(themodel)),
                                        AIC = as.numeric(AIC(themodel)),
                                        BIC = as.numeric(BIC(themodel)),
                                        MarR2 = as.numeric(r.squaredGLMM(themodel)[1]),
                                        ConR2 = as.numeric(r.squaredGLMM(themodel)[2]))
                output <- rbind(output, outputAdd)
                tTable <- data.frame(summary(themodel)$tTable)
                prevvari <- reducedIDV
                reducedIDV <- row.names(tTable)[tTable$p.value<0.05 & row.names(tTable) != "(Intercept)"]
                rm(tTable, outputAdd, reducedFomu, themodel)
              }
            }
            rm(i, prevvari)
            allsignificantV <- reducedIDV
            # drop one variable from all significant variables
            if(length(allsignificantV) > 1){
              tempV <- data.frame(t(combinat::combn(reducedIDV, (length(reducedIDV)-1))))
              for(i in 1:nrow(tempV)){
                reducedIDV <- as.character(tempV[i,])
                reducedFomu <- paste(DV, "~", reducedIDV[1], sep = "")
                if(length(reducedIDV)>=2){
                  for(j in 2:length(reducedIDV)){
                    reducedFomu <- paste(reducedFomu, "+", reducedIDV[j], sep = "")
                  }
                  rm(j)
                }
                themodelForm <- paste("themodel <- lme(", reducedFomu, ",",
                                      randomArg, ",",
                                      "data = data,
                                      control = lmeControl(maxIter=10000, msMaxIter = 10000))", sep = "")
                eval(parse(text=themodelForm))
                outputModel[[paste("ReducedModel", nrow(output)+1, sep = "")]] <- themodel
                rm(themodelForm)
                outputAdd <- data.table(Model = paste("ReducedModel", nrow(output)+1, sep = ""),
                                        Formula = reducedFomu,
                                        DIC = as.numeric(DIC(themodel)),
                                        AIC = as.numeric(AIC(themodel)),
                                        BIC = as.numeric(BIC(themodel)),
                                        MarR2 = as.numeric(r.squaredGLMM(themodel)[1]),
                                        ConR2 = as.numeric(r.squaredGLMM(themodel)[2]))
                output <- rbind(output, outputAdd)
                rm(outputAdd, reducedIDV, reducedFomu)
              }
              rm(i)
            }
            # add one variable to all significant variables from non-significant variables
            nonsignificantV <- allIDV[!(allIDV %in% allsignificantV)]
            for(addV in nonsignificantV){
              reducedIDV <- c(allsignificantV, addV)
              reducedFomu <- paste(DV, "~", reducedIDV[1], sep = "")
              if(length(reducedIDV)>=2){
                for(j in 2:length(reducedIDV)){
                  reducedFomu <- paste(reducedFomu, "+", reducedIDV[j], sep = "")
                }
                rm(j)
              }
              
              themodelForm <- paste("themodel <- lme(", reducedFomu, ",",
                                    randomArg, ",",
                                    "data = data,
                                    control = lmeControl(maxIter=10000, msMaxIter = 10000))", sep = "")
              eval(parse(text=themodelForm))
              outputModel[[paste("ReducedModel", nrow(output)+1, sep = "")]] <- themodel
              rm(themodelForm)
              outputAdd <- data.table(Model = paste("ReducedModel", nrow(output)+1, sep = ""),
                                      Formula = reducedFomu,
                                      DIC = as.numeric(DIC(themodel)),
                                      AIC = as.numeric(AIC(themodel)),
                                      BIC = as.numeric(BIC(themodel)),
                                      MarR2 = as.numeric(r.squaredGLMM(themodel)[1]),
                                      ConR2 = as.numeric(r.squaredGLMM(themodel)[2]))
              output <- rbind(output, outputAdd)
              rm(outputAdd, reducedFomu, reducedIDV, themodel)
            }
            rm(allsignificantV, nonsignificantV)
            return(invisible(list(modelSummary = output, modelOutput = outputModel)))
          })


#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(data = "data.frame",
                                DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                randomMethod = "character",
                                randomTerm = "character",
                                slopeTerm = "character"),
          definition = function(data,
                                DV,
                                IDV,
                                maxInteraction,
                                randomMethod,
                                randomTerm,
                                slopeTerm){
            return(data = data.table(data), DV, IDV, maxInteraction, randomMethod,
                   randomTerm, slopeTerm)
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(data = "data.table",
                                DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                randomMethod = "character",
                                randomTerm = "character",
                                slopeTerm = "character"),
          definition = function(data,
                                DV,
                                IDV,
                                randomMethod,
                                randomTerm,
                                slopeTerm){
            return(data, DV, IDV, maxInteraction = 1, randomMethod,
                   randomTerm, slopeTerm)
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(data = "data.table",
                                DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                randomMethod = "missing",
                                randomTerm = "character",
                                slopeTerm = "character"),
          definition = function(data,
                                DV,
                                IDV,
                                maxInteraction,
                                randomTerm,
                                slopeTerm){
            return(data, DV, IDV, maxInteraction, randomMethod = "intercept",
                   randomTerm, slopeTerm)
          })



#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(data = "data.frame",
                                DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                randomMethod = "character",
                                randomTerm = "character",
                                slopeTerm = "character"),
          definition = function(data,
                                DV,
                                IDV,
                                randomMethod,
                                randomTerm,
                                slopeTerm){
            return(data = data.table(data), DV, IDV, maxInteraction = 1, 
                   randomMethod,
                   randomTerm, slopeTerm)
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(data = "data.frame",
                                DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                randomMethod = "missing",
                                randomTerm = "character",
                                slopeTerm = "character"),
          definition = function(data,
                                DV,
                                IDV,
                                maxInteraction,
                                randomTerm,
                                slopeTerm){
            return(data = data.table(data), DV, IDV, maxInteraction,
                   randomMethod = "intercept",
                   randomTerm, slopeTerm)
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(data = "data.table",
                                DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                randomMethod = "missing",
                                randomTerm = "character",
                                slopeTerm = "character"),
          definition = function(data,
                                DV,
                                IDV,
                                randomTerm,
                                slopeTerm){
            return(data = data.table(data), DV, IDV, maxInteraction = 1,
                   randomMethod = "intercept",
                   randomTerm, slopeTerm)
          })

