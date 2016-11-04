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
                    distanceWeight) {
             standardGeneric("HeghyiCICalculation")
           })

#' @export
#' @rdname HeghyiCICalculation
setMethod(
  "HeghyiCICalculation",
  signature = c(data = "data.table",
                maxRadius = "numeric",
                sizeIndex = "character",
                distanceWeight = "numeric"),
  definition = function(data,
                        maxRadius,
                        sizeIndex,
                        distanceWeight){
    # calcuate coordination of each tree
    data[, ':='(coordX = sin(Angle*pi/180)*Distance,
                coordY = cos(Angle*pi/180)*Distance)]
    years <- sort(unique(data$Year))
    output <- data.table(PlotNumber = character(), TreeNumber = character(),
                         Year = numeric(), H = numeric(), IntraH = numeric(), InterH = numeric())
    
    for(indiyear in years){
      yeardata <- data[Year == indiyear,]
      yeardata[,temptreeno:=1:length(coordX), by = PlotNumber]
      for(i in 1:max(yeardata$temptreeno)){
        if(sizeIndex == "DBH"){
          targettrees <- yeardata[temptreeno == i, .(PlotNumber, TreeNumber,  temptreeno,
                                                     toSpecies = Species,toX = coordX,
                                                     toY = coordY, Distance, FocalSize = DBH)] 
          surroundingTrees <- yeardata[temptreeno != i, .(PlotNumber, NeighborSize = DBH, coordX, coordY, Species)]
        } else if(sizeIndex == "Biomass"){
          targettrees <- yeardata[temptreeno == i, .(PlotNumber, TreeNumber,  temptreeno,
                                                     toSpecies = Species,toX = coordX,
                                                     toY = coordY, Distance, FocalSize = Biomass)] 
          surroundingTrees <- yeardata[temptreeno != i, .(PlotNumber, NeighborSize = Biomass, coordX, coordY, Species)]
        } else {
          stop("Please specify sizeIndex from one of DBH or Biomass.")
        }
        surroundingTrees <- setkey(surroundingTrees, PlotNumber)[setkey(targettrees, PlotNumber),
                                                                 nomatch = 0]
        
        
        surroundingTrees[,XYDistance := (((coordX-toX)^2+(coordY-toY)^2)^0.5+0.1)]  
        surroundingTrees <- surroundingTrees[XYDistance <= maxRadius,]
        surroundingTrees_IntraSpecies <- surroundingTrees[Species == toSpecies,]
        
        totalHtable <- surroundingTrees[
          , .(tempH = sum(NeighborSize/(FocalSize*(XYDistance^distanceWeight)))),
          by = PlotNumber]
        if(nrow(surroundingTrees_IntraSpecies) > 0){
          IntraHtable <- surroundingTrees_IntraSpecies[
            , .(tempIntraH = sum(NeighborSize/(FocalSize*(XYDistance^distanceWeight)))),
            by = PlotNumber]
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
                        totalHtable[, .(PlotNumber, TreeNumber, Year = indiyear, 
                                        H = tempH*500/overLapArea,
                                        IntraH = tempIntraH*500/overLapArea,
                                        InterH = tempInterH*500/overLapArea)])
      }
      cat(indiyear, "is done. \n")
    }
    return(output)  
  })
