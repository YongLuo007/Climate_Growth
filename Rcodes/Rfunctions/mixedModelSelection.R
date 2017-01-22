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
#' @param ICCut numeric specify the cut point                         
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
#' @importFrom parallel makeCluster clusterExport
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
                    ICCut,
                    ...) {
             standardGeneric("mixedModelSelection")
           })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                ICTerm = "character",
                                ICCut = "numeric"),
          definition = function(DV,
                                IDV,
                                maxInteraction,
                                ICTerm,
                                ICCut,
                                ...){
            allIDV <- IDV
            output <- data.table(From = character(),
                                 Execution = character(),
                                 To = character(),
                                 NofIDV = numeric(),
                                 IC = numeric(),
                                 deltaIC = numeric())
            if(maxInteraction > 1){
              for(i in 2:maxInteraction){
                tempV <- data.frame(t(combinat::combn(IDV, i)), stringsAsFactors = FALSE)
                interactions <- as.character(tempV[,"X1"])
                for(j in 2:i){
                  interactions <- c(paste(interactions, ":",
                                          as.character(tempV[,paste("X", j, sep = "")]),
                                          sep = ""))
                }
                allIDV <- c(allIDV, interactions)
              }
            }
            # full model
            modelFormula <- paste(DV, "~", paste(allIDV, collapse = "+"), sep = "")
            themodel <- nlme::lme(as.formula(modelFormula),...)
            tTable <- data.frame(summary(themodel)$tTable)
            reducedIDV <- row.names(tTable)[tTable$p.value < 0.05 & row.names(tTable) != "(Intercept)"]
            outputAdd <- data.table(From = "NONE",
                                    Execution = "NONE",
                                    To = "Full",
                                    NofIDV = length(allIDV),
                                    IC = getIC(model = themodel, x = ICTerm),
                                    deltaIC = NA_real_)
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
                  BaseIC <- getIC(themodel, x = ICTerm)
                  outputAdd <- data.table(From = "Full",
                                          Execution = "All significant IDV",
                                          To = "Base",
                                          NofIDV = length(reducedIDV),
                                          IC = BaseIC,
                                          deltaIC = BaseIC-output$IC[1])
                  output <- rbind(output, outputAdd)
                  rm(outputAdd)
                }
                rm(tTable)
              }
            }
            rm(i, prevvari, reducedFomu)
            BaseIDV <- reducedIDV
            BaseModel <- themodel
            currentIDV <- BaseIDV
            
            previousModelName <- "Base"
            previousIC <- getIC(BaseModel, ICTerm)
            cl <- parallel::makeCluster(parallel::detectCores()-1)
            parallel::clusterExport(cl, c("getICformFomula", "lme", "AIC", "getIC"))
            reduceMark <- 1
            expendMark <- 1
            
            previousIDV <- character()
            for(i in 1:length(allIDV)){
              if(length(previousIDV) != length(currentIDV)){
                previousIDV <- currentIDV
                reducedFomus <- c()
                reducedIDV <- c()
                for(j in 1:length(currentIDV)){
                  reducedFomus <- c(reducedFomus, paste("logY~", paste(currentIDV[-j], collapse = "+"), sep = ""))
                  reducedIDV <- c(reducedIDV, currentIDV[j])
                }
                outputAdd_reduced <- data.table(From = previousModelName, 
                                                Execution = paste("-", reducedIDV, sep = ""),
                                                To = paste("Reduced-", reduceMark, "-", 1:length(currentIDV), sep = ""),
                                                NofIDV = length(currentIDV)-1,                 
                                                IC = unlist(parLapply(cl, reducedFomus, 
                                                                      function(y) getICformFomula(fomula = y, ICTerm, ...))))
                reduceMark <- reduceMark+1
                outputAdd_reduced[,':='(deltaIC = IC - previousIC, modelindex = 1:length(currentIDV))]
                bestModel <- outputAdd_reduced[deltaIC < ICCut & deltaIC == min(deltaIC),]
                if(nrow(bestModel) >= 1){
                  currentIDV <- previousIDV[-(bestModel[1,]$modelindex)]
                  previousIC <- bestModel[1,]$IC
                  previousModelName <- bestModel[1,]$To
                }
                output <- rbind(output, outputAdd_reduced[,modelindex:=NULL])
              }
            }  
            previousIDV <- character()
            for(i in 1:length(allIDV)){
              if(length(currentIDV) != length(previousIDV)){
                previousIDV <- currentIDV
                currentIDV <- unlist(lapply(lapply(lapply(currentIDV,  function(x){unlist(strsplit(x, ":"))}), sort),
                                            function(x){paste(x, collapse = ":")}))
                allIDV <- unlist(lapply(lapply(lapply(allIDV,  function(x){unlist(strsplit(x, ":"))}), sort),
                                        function(x){paste(x, collapse = ":")}))
                addedIDVs <- allIDV[!(allIDV %in% currentIDV)]
                if(length(addedIDVs)>0){
                  expandedFomus <- c()
                  expandedIDV <- c()
                  
                  for(j in 1:length(addedIDVs)){
                    expandedFomus <- c(expandedFomus, paste("logY~", paste(c(currentIDV, addedIDVs[j]), collapse = "+"), sep = ""))
                    expandedIDV <- c(expandedIDV, addedIDVs[j])
                  }
                  outputAdd_expanded <- data.table(From = previousModelName, 
                                                   Execution = paste("+", expandedIDV, sep = ""),
                                                   To = paste("Expanded-", expendMark, "-", 1:length(addedIDVs), sep = ""),
                                                   NofIDV = length(currentIDV)+1, 
                                                   IC = unlist(parLapply(cl, expandedFomus, 
                                                                         function(y) getICformFomula(fomula = y, ICTerm, ...))))
                  expendMark <- expendMark+1
                  outputAdd_expanded[,':='(deltaIC = IC - previousIC, modelindex = 1:length(addedIDVs))]
                  bestModel <- outputAdd_expanded[deltaIC < -ICCut & deltaIC == min(deltaIC),]
                  if(nrow(bestModel) >= 1){
                    currentIDV <- c(previousIDV, addedIDVs[(bestModel[1,]$modelindex)])
                    previousIC <- bestModel[1,]$IC
                    previousModelName <- bestModel[1,]$To
                  }
                  output <- rbind(output, outputAdd_expanded[,modelindex:=NULL])
                }

              }
            }
            
            output[, bestModel := "No"]
            output[To == previousModelName, bestModel := "Yes"]
            stopCluster(cl)
            return(invisible(list(modelSummary = output, 
                                  bestIDV = currentIDV)))
          })

getICformFomula <- function(fomula, ICTerm, ...) {
  themodel <- nlme::lme(fixed = as.formula(fomula), ...)
  return(getIC(themodel, x = ICTerm))
}

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




#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                ICTerm = "character",
                                ICCut = "numeric"),
          definition = function(DV,
                                IDV,
                                ICTerm,
                                ICCut,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction = 1, ICTerm, ICCut, ...))
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                ICTerm = "missing",
                                ICCut = "numeric"),
          definition = function(DV,
                                IDV,
                                maxInteraction,
                                ICCut,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction, ICTerm = "AIC", ICCut, ...))
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                ICTerm = "character",
                                ICCut = "missing"),
          definition = function(DV,
                                IDV,
                                maxInteraction,
                                ICTerm,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction, ICTerm, ICCut = 2, ...))
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                ICTerm = "missing",
                                ICCut = "numeric"),
          definition = function(DV,
                                IDV,
                                ICCut,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction = 1, ICTerm = "AIC", ICCut, ...))
          })


#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                ICTerm = "character",
                                ICCut = "missing"),
          definition = function(DV,
                                IDV,
                                ICTerm,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction = 1, ICTerm, ICCut = 2, ...))
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "numeric",
                                ICTerm = "missing",
                                ICCut = "missing"),
          definition = function(DV,
                                IDV,
                                maxInteraction,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction, ICTerm = "AIC", ICCut = 2, ...))
          })

#' @export
#' @rdname mixedModelSelection
setMethod("mixedModelSelection",
          signature = signature(DV = "character",
                                IDV = "character",
                                maxInteraction = "missing",
                                ICTerm = "missing",
                                ICCut = "missing"),
          definition = function(DV,
                                IDV,
                                ...){
            return(mixedModelSelection(DV, IDV, maxInteraction = 1, ICTerm = "AIC", ICCut = 2, ...))
          })
