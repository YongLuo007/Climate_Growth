#' the function to calculate intraspecific and interspecific hegyi competition index using DBH or biomass
#'
#' @param data data.table which must have PlotNumber TreeNumber and Year at which the trees were
#'                        measured. DBH must be present if using DBH, Biomass must be present if
#'                        using biomass to calculate competition index, Distance and Angle, and species
#' 
#' @param maxRadius numeric, the competition index will been calculated within this radius
#' 
#' @param sizeIndex character, choose DBH or Biomass to calculate competition
#' 
#'
#' @return a data table that has five columns, plotNumber, treeNumber, Year, IntraH and InterH
#' 
#' @importFrom data.table data.table ':='
#' @importFrom dplyr left_join '%>%' 
#'
#' @note no note
#'
#' @seealso no
#'
#' @export
#' @docType methods
#' @rdname HeghyiCICalculation
#'
#' @author Yong Luo
#'
setGeneric("HeghyiCICalculation",
           function(data,
                    maxRadius,
                    sizeIndex,
                    distanceWeight,
                    sizeWeight) {
             standardGeneric("HeghyiCICalculation")
           })

#' @export
#' @rdname HeghyiCICalculation
setMethod(
  "HeghyiCICalculation",
  signature = c(data = "data.table",
                maxRadius = "numeric",
                sizeIndex = "character",
                distanceWeight = "numeric",
                sizeWeight = "numeric"),
  definition = function(data,
                        maxRadius,
                        sizeIndex,
                        distanceWeight,
                        sizeWeight){
    # calcuate coordination of each tree
    alldata <- list()
    data[, ':='(coordX = sin(Angle*pi/180)*Distance,
                coordY = cos(Angle*pi/180)*Distance)]
    years <- sort(unique(data$Year))
    output <- data.table(PlotNumber = character(), TreeNumber = character(),
                         Year = numeric(), H = numeric(), IntraH = numeric(), InterH = numeric())
    for(indiyear in years){
      yeardata <- data[Year == indiyear,]
      yeardata[,temptreeno:=1:length(coordX), by = PlotNumber]
      alldata[[paste(indiyear)]] <- yeardata
    }
    newEnv <- new.env()
    newEnv$alldata <- alldata
    newEnv$sizeIndex <- sizeIndex
    newEnv$maxRadius <- maxRadius
    newEnv$distanceWeight <- distanceWeight
    newEnv$sizeWeight <- sizeWeight
    rm(data)
    cl <- parallel::makeCluster(detectCores()-1)
    parallel::clusterExport(cl, c("alldata", "sizeIndex", "maxRadius", 
                                  "mainFunc", "distanceWeight", "sizeWeight"), envir = newEnv)
    parallel::clusterExport(cl, c("data.table", "setkey", "%>%"))
    alloutput <- parLapply(cl, alldata, function(x) mainFunc(yeardata = x, sizeIndex = sizeIndex,
                                                             maxRadius = maxRadius, 
                                                             distanceWeight = distanceWeight,
                                                             sizeWeight = sizeWeight))
    stopCluster(cl)
    for(i in 1:length(alloutput)){
      output <- rbind(output, alloutput[[i]])
    }
    return(output)  
  })


mainFunc <- function(yeardata, sizeIndex, maxRadius, distanceWeight, sizeWeight){
  output <- data.table(PlotNumber = character(), TreeNumber = character(),
                       Year = numeric(), H = numeric(), IntraH = numeric(), InterH = numeric())
  for(i in 1:max(yeardata$temptreeno)){
    if(sizeIndex == "DBH"){
      sizeRangedata <- yeardata[,.(minSize = min(DBH), maxSize = max(DBH)), by = PlotNumber]
      targettrees <- yeardata[temptreeno == i, .(PlotNumber, TreeNumber,  temptreeno,
                                                 toSpecies = Species,toX = coordX,
                                                 toY = coordY, Distance, FocalSize = DBH)] 
      surroundingTrees <- yeardata[temptreeno != i, .(PlotNumber, NeighborSize = DBH, coordX, coordY, Species)]
    } else if(sizeIndex == "Biomass"){
      sizeRangedata <- yeardata[,.(minSize = min(Biomass), maxSize = max(Biomass)), by = PlotNumber]
      targettrees <- yeardata[temptreeno == i, .(PlotNumber, TreeNumber,  temptreeno,
                                                 toSpecies = Species,toX = coordX,
                                                 toY = coordY, Distance, FocalSize = Biomass)] 
      surroundingTrees <- yeardata[temptreeno != i, .(PlotNumber, NeighborSize = Biomass, coordX, coordY, Species)]
    } else {
      stop("Please specify sizeIndex from one of DBH or Biomass.")
    }
    surroundingTrees <- setkey(surroundingTrees, PlotNumber)[setkey(targettrees, PlotNumber),
                                                             nomatch = 0]
    surroundingTrees <- setkey(surroundingTrees, PlotNumber)[setkey(sizeRangedata, PlotNumber),
                                                             nomatch = 0]
    
    surroundingTrees[,XYDistance := (((coordX-toX)^2+(coordY-toY)^2)^0.5+0.1)]  
    surroundingTrees <- surroundingTrees[XYDistance <= maxRadius,]
    surroundingTrees_IntraSpecies <- surroundingTrees[Species == toSpecies,]
    
    totalHtable <- surroundingTrees[,.(tempH = ((exp((maxSize-FocalSize)/(maxSize-minSize)))^sizeWeight)/FocalSize*
                                         sum(NeighborSize/(XYDistance^distanceWeight))),
                                    by = PlotNumber] %>%
      unique(., by = "PlotNumber")
    if(nrow(surroundingTrees_IntraSpecies) > 0){
      IntraHtable <- surroundingTrees_IntraSpecies[,.(tempIntraH = ((exp((maxSize-FocalSize)/(maxSize-minSize)))^sizeWeight)/FocalSize*
                                                        sum(NeighborSize/(XYDistance^distanceWeight))),
                                                   by = PlotNumber] %>%
        unique(., by = "PlotNumber")
      totalHtable <- dplyr::left_join(totalHtable, IntraHtable, by = "PlotNumber") %>%
        data.table
      totalHtable[is.na(tempIntraH), tempIntraH := 0]
    } else {
      totalHtable[,tempIntraH := 0]
    }
    totalHtable[, tempInterH := tempH - tempIntraH]
    totalHtable <- setkey(totalHtable, PlotNumber)[setkey(targettrees[,.(PlotNumber, TreeNumber, 
                                                                         Distance)],
                                                          PlotNumber), nomatch = 0]
    totalHtable[, overLapArea := 2*(maxRadius^2)*acos(0.5*Distance/maxRadius)-0.5*Distance*((4*(maxRadius^2)-Distance^2)^0.5)]
    
    output <- rbind(output, 
                    totalHtable[, .(PlotNumber, TreeNumber, Year = yeardata$Year[1], 
                                    H = tempH*500/overLapArea,
                                    IntraH = tempIntraH*500/overLapArea,
                                    InterH = tempInterH*500/overLapArea)])
  }
  return(output)
}



