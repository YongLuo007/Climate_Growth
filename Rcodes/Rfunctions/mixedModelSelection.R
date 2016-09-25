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
                    ...) {
             standardGeneric("mixedModelSelection")
           })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric"),
          definition = function(DV,
                                IDV,
                                maxInteraction,
                                ...){
            allIDV <- IDV
            output <- data.table(Model = character(),
                                 IDV_Base = character(),
                                 Direction = character(),
                                 IDV_Processed = character(),
                                 IDV_Length = numeric(),
                                 All_Significant = logical(),
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
            themodel <- nlme::lme(as.formula(modelFormula),...)
            outputModel[["Full"]] <- themodel
            tTable <- data.frame(summary(themodel)$tTable)
            reducedIDV <- row.names(tTable)[tTable$p.value < 0.05 & row.names(tTable) != "(Intercept)"]
            outputAdd <- data.table(Model = "Full",
                                    IDV_Base = "NONE",
                                    Direction = "NONE",
                                    IDV_Processed = "NONE",
                                    IDV_Length = length(allIDV),
                                    All_Significant = length(allIDV)==length(reducedIDV),
                                    DIC = as.numeric(MuMIn::DIC(themodel)),
                                    AIC = as.numeric(AIC(themodel)),
                                    BIC = as.numeric(BIC(themodel)),
                                    MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                                    ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
            output <- rbind(output, outputAdd)
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
                prevvari <- reducedIDV
                themodel <- nlme::lme(as.formula(reducedFomu),...)
                tTable <- data.frame(summary(themodel)$tTable)
                reducedIDV <- row.names(tTable)[tTable$p.value < 0.05 & row.names(tTable) != "(Intercept)"]
                if(length(reducedIDV) == length(prevvari)){
                  outputModel[["BaseModel"]] <- themodel
                  BaseIDV <- paste(reducedIDV, collapse = ", ")
                  outputAdd <- data.table(Model = "BaseModel",
                                          IDV_Base = BaseIDV,
                                          Direction = "BaseModel",
                                          IDV_Processed = "NONE",
                                          IDV_Length = length(reducedIDV),
                                          All_Significant = TRUE,
                                          DIC = as.numeric(MuMIn::DIC(themodel)),
                                          AIC = as.numeric(AIC(themodel)),
                                          BIC = as.numeric(BIC(themodel)),
                                          MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                                          ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
                  output <- rbind(output, outputAdd)
                  rm(outputAdd)
                }
                rm(tTable, reducedFomu, themodel)
              }
            }
            rm(i, prevvari)
            
            allsignificantV <- reducedIDV
            # drop one variable from all significant variables
            if(length(allsignificantV) > 1){
              for(k in 1:(length(allsignificantV)-1)){
                tempV <- data.frame(t(combinat::combn(allsignificantV,
                                                      (length(allsignificantV)-k))),
                                    stringsAsFactors = FALSE)
                for(i in 1:nrow(tempV)){
                  reducedIDV <- as.character(tempV[i,])
                  reducedFomu <- paste(DV, "~", reducedIDV[1], sep = "")
                  if(length(reducedIDV)>=2){
                    for(j in 2:length(reducedIDV)){
                      reducedFomu <- paste(reducedFomu, "+", reducedIDV[j], sep = "")
                    }
                    rm(j)
                  }
                  themodel <- nlme::lme(as.formula(reducedFomu),...)
                  outputModel[[paste("ReducedModel", nrow(output)-1, sep = "")]] <- themodel
                  droppedV <- paste(allsignificantV[!(allsignificantV %in% reducedIDV)], collapse=", ")
                  tTable <- data.frame(summary(themodel)$tTable)
                  sigV <- row.names(tTable)[tTable$p.value < 0.05 & row.names(tTable) != "(Intercept)"]
                  outputAdd <- data.table(Model = paste("ReducedModel", nrow(output)-1, sep = ""),
                                          IDV_Base = BaseIDV,
                                          Direction = "drop",
                                          IDV_Processed = paste(droppedV, collapse = ", "),
                                          IDV_Length = length(reducedIDV),
                                          All_Significant = length(reducedIDV) == length(sigV),
                                          DIC = as.numeric(MuMIn::DIC(themodel)),
                                          AIC = as.numeric(AIC(themodel)),
                                          BIC = as.numeric(BIC(themodel)),
                                          MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                                          ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
                  output <- rbind(output, outputAdd)
                  rm(outputAdd, reducedFomu, reducedIDV, tTable)
                }
                rm(i, tempV)
              }
            }
            
            # add one variable to all significant variables from non-significant variables
            nonsignificantV <- allIDV[!(allIDV %in% allsignificantV)]
            nonsignificantVLength <- length(nonsignificantV)
            k <- 1
            significantList <- list()
            AICS <- c()
            for(addV in nonsignificantV){
              reducedIDV <- c(allsignificantV, addV)
              reducedFomu <- paste(DV, "~", reducedIDV[1], sep = "")
              if(length(reducedIDV)>=2){
                for(j in 2:length(reducedIDV)){
                  reducedFomu <- paste(reducedFomu, "+", reducedIDV[j], sep = "")
                }
                rm(j)
              }
              themodel <- nlme::lme(as.formula(reducedFomu),...)
              outputModel[[paste("ExpandedModel", nrow(output)-1, sep = "")]] <- themodel
              tTable <- data.frame(summary(themodel)$tTable)
              sigV <- row.names(tTable)[tTable$p.value < 0.05 & row.names(tTable) != "(Intercept)"]
              outputAdd <- data.table(Model = paste("ExpandedModel", nrow(output)-1, sep = ""),
                                      IDV_Base = BaseIDV,
                                      Direction = "add",
                                      IDV_Processed = paste(addV, collapse = ", "),
                                      IDV_Length = length(reducedIDV),
                                      All_Significant = length(reducedIDV) == length(sigV),
                                      DIC = as.numeric(MuMIn::DIC(themodel)),
                                      AIC = as.numeric(AIC(themodel)),
                                      BIC = as.numeric(BIC(themodel)),
                                      MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                                      ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
              if(length(reducedIDV) == length(sigV)){
                significantList[[k]] <- reducedIDV
                k <- k+1
                AICS <- c(AICS, as.numeric(AIC(themodel)))
              }
              output <- rbind(output, outputAdd)
              rm(outputAdd, reducedFomu, reducedIDV, themodel)
            }
            rm(allsignificantV, nonsignificantV)
            for(i in 1:(nonsignificantVLength-1)){
              if(k > 1){
                nonsignificantV <- allIDV[!(allIDV %in% significantList[[which.min(AICS)]])]
                allsignificantV <- significantList[[which.min(AICS)]]
                addedIDVs <- allsignificantV[!(allsignificantV %in% unlist(strsplit(BaseIDV, ", ", fixed = TRUE)))]
                k <- 1
                significantList <- list()
                for(addV in nonsignificantV){
                  reducedIDV <- c(allsignificantV, addV)
                  reducedFomu <- paste(DV, "~", reducedIDV[1], sep = "")
                  if(length(reducedIDV)>=2){
                    for(j in 2:length(reducedIDV)){
                      reducedFomu <- paste(reducedFomu, "+", reducedIDV[j], sep = "")
                    }
                    rm(j)
                  }
                  themodel <- nlme::lme(as.formula(reducedFomu),...)
                  outputModel[[paste("ExpandedModel", nrow(output)-1, sep = "")]] <- themodel
                  tTable <- data.frame(summary(themodel)$tTable)
                  sigV <- row.names(tTable)[tTable$p.value < 0.05 & row.names(tTable) != "(Intercept)"]
                  outputAdd <- data.table(Model = paste("ExpandedModel", nrow(output)-1, sep = ""),
                                          IDV_Base = BaseIDV,
                                          Direction = "add",
                                          IDV_Processed = paste(c(addedIDVs, addV), collapse = ", "),
                                          IDV_Length = length(reducedIDV),
                                          All_Significant = length(reducedIDV) == length(sigV),
                                          DIC = as.numeric(MuMIn::DIC(themodel)),
                                          AIC = as.numeric(AIC(themodel)),
                                          BIC = as.numeric(BIC(themodel)),
                                          MarR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[1]),
                                          ConR2 = as.numeric(MuMIn::r.squaredGLMM(themodel)[2]))
                  if(length(reducedIDV) == length(sigV)){
                    significantList[[k]] <- reducedIDV
                    k <- k+1
                  }
                  output <- rbind(output, outputAdd)
                  rm(outputAdd, reducedFomu, reducedIDV, themodel)
                }
                rm(allsignificantV, nonsignificantV)
              }
            }
            return(invisible(list(modelSummary = output, modelOutput = outputModel)))
          })


#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "missing"),
          definition = function(DV,
                                IDV,
                                ...){
            return(DV, IDV, maxInteraction = 1, ...)
          })
