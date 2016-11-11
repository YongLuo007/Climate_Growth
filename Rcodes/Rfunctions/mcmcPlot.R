#' To produce the trace plot and density plot by a vairiable name
#'
#' @param x mcmc.list, a mcmc.list contains one or several chains
#' 
#' @param variableName character, specify the variable name that must be sampled.
#' 
#'
#' @return no object returned
#' 
#' @importFrom data.table data.table ':='
#' @importFrom dplyr left_join '%>%' 
#' @importFrom ggplot2 ggplot ggplotGrob
#' @importFrom coda nchain
#' @importFrom gridExtra grid.arrange
#'
#' @note no note
#'
#' @seealso no
#'
#' @example 
#' \dontrun{
#' library(coda)
#' #load 3 mcmc.lists with each list contains one chain
#' for(k in 1:3){
#' load(file.path("~", "GitHub", "Climate_Growth", "Results", "JAGS", "TA", 
#'                paste("chain", k, 
#'                      ".RData", sep = "")))
#' if(k == 1){
#'   allcoda <- get(paste("coda", k, sep = ""))
#' } else {
#'   allcoda[k] <- get(paste("coda", k, sep = ""))
#' }
#' }
#' }
#' a <- mcmcPlot(allcoda, "a")
#'
#' @export
#' @docType methods
#' @importClassesFrom coda mcmc.list
#' @rdname mcmcPlot
#'
#' @author Yong Luo
#'
setGeneric("mcmcPlot",
           function(x,
                    variableName) {
             standardGeneric("mcmcPlot")
           })

#' @export
#' @rdname mcmcPlot
setMethod(
  "mcmcPlot",
  signature = c(x = "mcmc.list",
                variableName = "character"),
  definition = function(x,
                        variableName){
    codatable <- data.table(as.matrix(x))[,':='(Time = time(x), 
                                                chain = sort(rep(1:nchain(x), length(time(x)))))]
    notInVariable <- variableName[!(variableName %in% names(codatable))]
    if(length(notInVariable)>0){
      cat("Variables", paste(notInVariable, collapse = ", "), "not in mcmc.list. The rests will be plotted.")
      variableName <- variableName[!(variableName %in% notInVariable)]
    }
    if(length(variableName) == 0){
      stop("No valid variable names provided, please verify your variable name.")
    }
    conNum <- which(names(codatable) %in% c(variableName))
    subcodas <- x[,conNum, drop = TRUE]
    if(nchain(x) > 1){
      gelmantest <- gelman.diag(subcodas, multivariate = FALSE)
      if(length(variableName) == 1){
        gelmanResults <- data.table(gelmantest$psrf)[
          ,.(variable = variableName, GelmenTest = as.character(sprintf("%0.2f", round(`Point est.`, digits = 2))))]
      } else {
        gelmanResults <- data.table(gelmantest$psrf, keep.rownames = TRUE)[
          , .(variable = rn, GelmenTest = as.character(sprintf("%0.2f", round(`Point est.`, digits = 2))))]
      }
    } else {
      gelmanResults <- data.table(variable = variableName, GelmenTest = "NA")
    }
   
    codatable <- codatable[ ,names(codatable) %in% c(variableName, "chain", "Time"), with = FALSE]
    codaLongForm <- melt(codatable, id.vars = c("chain", "Time"), measure.vars = variableName,
                         value.name = "Value")
    densityTable <- data.table(chain = numeric(), variable = character(), Density = numeric(), 
                               Value = numeric())
    for(j in variableName){
      for(i in 1:nchain(x)){
        tempdensity <- density(codaLongForm[chain == i & variable == j, ]$Value,
                               from = min(codaLongForm[chain == i & variable == j, ]$Value),
                               to = max(codaLongForm[chain == i & variable == j, ]$Value))
        densityTable <- rbind(densityTable,
                              data.table(chain = i, variable = j, Density = 0,
                                         Value = min(tempdensity$x)),
                              data.table(chain = i, variable = j,
                                         Density = tempdensity$y,
                                         Value = tempdensity$x),
                              data.table(chain = i, variable = j, Density = 0,
                                         Value = max(tempdensity$x)))
        rm(tempdensity)
      }
      tempdensity <- density(codaLongForm[variable == j, ]$Value,
                             from = min(codaLongForm[variable == j, ]$Value),
                             to = max(codaLongForm[variable == j, ]$Value))
      densityTable <- rbind(densityTable,
                            data.table(chain = 0, variable = j, Density = 0,
                                       Value = min(tempdensity$x)),
                            data.table(chain = 0, variable = j,
                                       Density = tempdensity$y,
                                       Value = tempdensity$x),
                            data.table(chain = 0, variable = j, Density = 0,
                                       Value = max(tempdensity$x)))
    }
    codaLongForm[,':='(variable = factor(variable, levels = variableName),
                       Chain = factor(chain, levels = 1:nchain(x),
                                      labels = paste("Chain ", 1:nchain(x), sep = "")))]
    densityTable[,':='(variable = factor(variable, levels = variableName),
                       Chain = factor(chain, levels = 1:nchain(x),
                                      labels = paste("Chain ", 1:nchain(x), sep = "")))]
    addtoGelman1 <- codaLongForm[,.(y = max(Value) + abs(max(Value)-min(Value))/15, 
                                    Time = min(Time)+(max(Time)-min(Time))/4),
                                 by = variable]
    addtoGelman2 <- densityTable[,.(Density = min(Density)+(max(Density)-min(Density))/2),
                                 by = variable]
    
    texts <- setkey(gelmanResults, variable)[setkey(addtoGelman1, variable), nomatch = 0]
    texts <- setkey(texts, variable)[setkey(addtoGelman2, variable), nomatch = 0]
    texts[,':='(text1 = paste("No. chains: ", nchain(x), "; Gelman test: ", GelmenTest, sep = ""),
                        text2 = "Black line is for all chains")]
    
    tracePlot <- ggplot(data = codaLongForm, aes(x = Time, y = Value))+
      geom_line(aes(group = Chain, col = Chain))+
      geom_text(data = texts, aes(x = Time, y = y, label = text1), hjust = 0)+
      facet_grid(variable~., scales = "free_y", switch = "y")+
      theme_bw()+
      theme(panel.grid = element_blank(),
            legend.position = "none",
            strip.text.y = element_text(size = rel(1.3)),
            strip.background = element_rect(fill = "white", colour = "white"),
            axis.title.y = element_blank())
    
    densityPlot <- ggplot(data = densityTable[chain != 0,], 
                          aes(x = Density, y = Value))+
      facet_grid(variable~., scales = "free_y")+
      geom_polygon(aes(group = Chain, col = Chain), 
                   fill = "white", alpha = 0)+
      geom_polygon(data = densityTable[chain == 0,], aes(x = Density, y = Value),
                   fill = "white", alpha = 0, col = "black", size = rel(1.1))+
      geom_vline(xintercept = 0, col = "white", size = rel(1.1))+
      geom_text(data = texts[variable == variableName[1], ], 
                aes(x = Density, y = y, label = text2))+
      theme_bw()+
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text = element_blank(),
            legend.position = "none")
    tracePlotTree <- ggplotGrob(tracePlot)
    densityPlotTree <- ggplotGrob(densityPlot)
    
    tracePlotTree$heights <- densityPlotTree$heights
    plotlayout <- cbind(c(1), c(1), c(2))
    allplots <- gridExtra::grid.arrange(tracePlotTree, densityPlotTree,
                                        layout_matrix = plotlayout)
    return(invisible(allplots))
  })
