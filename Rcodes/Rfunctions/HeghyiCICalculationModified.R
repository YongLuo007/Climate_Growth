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
                    sizeWeight,
                    assymetricScale,
                    testing) {
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
                sizeWeight = "numeric",
                assymetricScale = "character",
                testing = "logical"),
  definition = function(data,
                        maxRadius,
                        sizeIndex,
                        distanceWeight,
                        sizeWeight,
                        assymetricScale,
                        testing = FALSE){
    # calcuate coordination of each tree
    alldata <- list()
    data[, ':='(coordX = sin(Angle*pi/180)*Distance,
                coordY = cos(Angle*pi/180)*Distance)]
    years <- sort(unique(data$Year))
    for(indiyear in years){
      yeardata <- data[Year == indiyear,]
      yeardata[,temptreeno:=1:length(coordX), by = PlotNumber]
      alldata[[paste(indiyear)]] <- yeardata
    }
    # browser()
    # a <- mainFunc(yeardata = alldata[[1]], sizeIndex = sizeIndex,
    #               maxRadius = maxRadius, 
    #               distanceWeight = distanceWeight,
    #               sizeWeight = sizeWeight,
    #               assymetricScale = assymetricScale,
    #               testing = testing)
    newEnv <- new.env()
    newEnv$alldata <- alldata
    newEnv$sizeIndex <- sizeIndex
    newEnv$maxRadius <- maxRadius
    newEnv$distanceWeight <- distanceWeight
    newEnv$sizeWeight <- sizeWeight
    newEnv$assymetricScale <- assymetricScale
    newEnv$testing <- testing
    rm(data)
    cl <- parallel::makeCluster(detectCores()-1)
    parallel::clusterExport(cl, c("alldata", "sizeIndex", "maxRadius", 
                                  "mainFunc", "distanceWeight", "sizeWeight",
                                  "assymetricScale", "testing"), envir = newEnv)
    parallel::clusterExport(cl, c("data.table", "setkey", "%>%", "dcast", "setcolorder", "unique"))
    alloutput <- parLapply(cl, alldata, function(x) mainFunc(yeardata = x, sizeIndex = sizeIndex,
                                                             maxRadius = maxRadius, 
                                                             distanceWeight = distanceWeight,
                                                             sizeWeight = sizeWeight,
                                                             assymetricScale = assymetricScale,
                                                             testing = testing))
    stopCluster(cl)
    rm(newEnv)
    for(i in 1:length(alloutput)){
      if(i == 1){
        output <- alloutput[[1]]$output
        verifyTableoutput <- alloutput[[1]]$verifyOutput
      } else {
        output <- rbind(output, alloutput[[i]]$output)
        verifyTableoutput <- rbind(verifyTableoutput, alloutput[[i]]$verifyOutput)
      }
    }
    if(!(testing)){
      return(list(output = output))
    } else {
      return(list(output = output, verifyOutput = verifyTableoutput))
    }  
  })


mainFunc <- function(yeardata, sizeIndex, maxRadius, distanceWeight,
                     sizeWeight, assymetricScale, testing){
  output <- data.table(PlotNumber = character(), TreeNumber = character(),
                       Year = numeric(), H = numeric(), IntraH = numeric(), InterH = numeric())
  weightTable <- data.table(expand.grid(distanceWeight = distanceWeight, sizeWeight = sizeWeight,
                                        stringsAsFactors = FALSE))
  weightTable[,competitionName := paste("DW", distanceWeight, "_SW", sizeWeight, sep = "")]
  for(i in 1:max(yeardata$temptreeno)){
    if(sizeIndex == "DBH"){
      sizeRangedata <- yeardata[,.(minSize = min(DBH), meanSize = mean(DBH), maxSize = max(DBH)), by = PlotNumber]
      targettrees <- yeardata[temptreeno == i, .(PlotNumber, TreeNumber,  temptreeno,
                                                 toSpecies = Species,
                                                 toAngle = Angle,
                                                 toDistance = Distance, FocalSize = DBH)] 
      surroundingTrees <- yeardata[temptreeno != i, .(PlotNumber, NeighborSize = DBH, Distance, Angle, Species)]
    } else if(sizeIndex == "Biomass"){
      sizeRangedata <- yeardata[,.(minSize = min(Biomass), meanSize = mean(Biomass), maxSize = max(Biomass)), by = PlotNumber]
      targettrees <- yeardata[temptreeno == i, .(PlotNumber, TreeNumber,  temptreeno,
                                                 toSpecies = Species,toAngle = Angle,
                                                 toDistance = Distance, FocalSize = Biomass)] 
      surroundingTrees <- yeardata[temptreeno != i, .(PlotNumber, neigborN = temptreeno,NeighborSize = Biomass, Distance, Angle, Species)]
    } else {
      stop("Please specify sizeIndex from one of DBH or Biomass.")
    }
    surroundingTrees <- setkey(surroundingTrees, PlotNumber)[setkey(targettrees, PlotNumber),
                                                             nomatch = 0]
    surroundingTrees[,':='(Angle = 180-toAngle+Angle)]
    surroundingTrees[,':='(toAngle = 180)]
    surroundingTrees[, ':='(coordX = sin(Angle*pi/180)*Distance,
                coordY = cos(Angle*pi/180)*Distance)]
    surroundingTrees[, ':='(toX = sin(toAngle*pi/180)*toDistance,
                            toY = cos(toAngle*pi/180)*toDistance)]
    surroundingTrees[, coordY:=coordY+abs(toY)]
    surroundingTrees[, toY := 0]
    
    surroundingTrees <- setkey(surroundingTrees, PlotNumber)[setkey(sizeRangedata, PlotNumber),
                                                             nomatch = 0]
    surroundingTrees[,XYDistance := (((coordX-toX)^2+(coordY-toY)^2)^0.5+0.1)]  
    surroundingTrees <- surroundingTrees[XYDistance <= maxRadius,]
    surroundingTrees <- setkey(surroundingTrees[,k:=1], k)[setkey(weightTable[,k:=1], k),
                                                           nomatch = NA, allow.cartesian = TRUE]
    surroundingTrees13 <- surroundingTrees[coordX/coordY <= tan(60*pi/180) &
                                              coordY >= 0 & tan(300*pi/180) <= coordX/coordY,]
    surroundingTrees_IntraSpecies <- surroundingTrees[Species == toSpecies,]
    surroundingTrees13_IntraSpecies <- surroundingTrees13[Species == toSpecies,]
    if(assymetricScale == "Rescale"){
    totalHtable <- surroundingTrees[,.(tempH = ((exp((maxSize-FocalSize)/(maxSize-minSize)))^sizeWeight)/FocalSize*
                                         sum(NeighborSize/(XYDistance^distanceWeight))),
                                    by = c("PlotNumber", "competitionName")] %>%
      unique(., by = c("PlotNumber", "competitionName"))
    totalHtable13 <- surroundingTrees13[,.(temp13H = ((exp((maxSize-FocalSize)/(maxSize-minSize)))^sizeWeight)/FocalSize*
                                            sum(NeighborSize/(XYDistance^distanceWeight))),
                                       by = c("PlotNumber", "competitionName")] %>%
      unique(., by = c("PlotNumber", "competitionName"))
    totalHtable <- dplyr::left_join(totalHtable, totalHtable13, by = c("PlotNumber", "competitionName")) %>%
      data.table
    totalHtable[is.na(temp13H), temp13H := 0]
    } else if (assymetricScale == "Relative"){
      totalHtable <- surroundingTrees[,.(tempH = ((meanSize/FocalSize)^sizeWeight)/FocalSize*
                                           sum(NeighborSize/(XYDistance^distanceWeight))),
                                      by = c("PlotNumber", "competitionName")] %>%
        unique(., by = c("PlotNumber", "competitionName"))
      totalHtable13 <- surroundingTrees13[,.(tempH13 = ((meanSize/FocalSize)^sizeWeight)/FocalSize*
                                           sum(NeighborSize/(XYDistance^distanceWeight))),
                                      by = c("PlotNumber", "competitionName")] %>%
        unique(., by = c("PlotNumber", "competitionName"))
      totalHtable <- dplyr::left_join(totalHtable, totalHtable13, by = c("PlotNumber", "competitionName")) %>%
        data.table
      totalHtable[is.na(temp13H), temp13H := 0]
    }
    if(nrow(surroundingTrees_IntraSpecies) > 0){
      if(assymetricScale == "Rescale"){
        IntraHtable <- surroundingTrees_IntraSpecies[,.(tempIntraH = ((exp((maxSize-FocalSize)/(maxSize-minSize)))^sizeWeight)/FocalSize*
                                             sum(NeighborSize/(XYDistance^distanceWeight))),
                                        by = c("PlotNumber", "competitionName")] %>%
          unique(., by = c("PlotNumber", "competitionName"))
      } else if (assymetricScale == "Relative"){
        IntraHtable <- surroundingTrees_IntraSpecies[,.(tempIntraH = ((meanSize/FocalSize)^sizeWeight)/FocalSize*
                                             sum(NeighborSize/(XYDistance^distanceWeight))),
                                        by = c("PlotNumber", "competitionName")] %>%
          unique(., by = c("PlotNumber", "competitionName"))
      }
      totalHtable <- dplyr::left_join(totalHtable, IntraHtable, by = c("PlotNumber", "competitionName")) %>%
        data.table
      totalHtable[is.na(tempIntraH), tempIntraH := 0]
    } else {
      totalHtable[,tempIntraH := 0]
    }
    if(nrow(surroundingTrees13_IntraSpecies) > 0){
      if(assymetricScale == "Rescale"){
        IntraHtable13 <- surroundingTrees13_IntraSpecies[,.(tempIntraH13 = ((exp((maxSize-FocalSize)/(maxSize-minSize)))^sizeWeight)/FocalSize*
                                                          sum(NeighborSize/(XYDistance^distanceWeight))),
                                                     by = c("PlotNumber", "competitionName")] %>%
          unique(., by = c("PlotNumber", "competitionName"))
      } else if (assymetricScale == "Relative"){
        IntraHtable13 <- surroundingTrees13_IntraSpecies[,.(tempIntraH13 = ((meanSize/FocalSize)^sizeWeight)/FocalSize*
                                                          sum(NeighborSize/(XYDistance^distanceWeight))),
                                                     by = c("PlotNumber", "competitionName")] %>%
          unique(., by = c("PlotNumber", "competitionName"))
      }
      totalHtable <- dplyr::left_join(totalHtable, IntraHtable13, by = c("PlotNumber", "competitionName")) %>%
        data.table
      totalHtable[is.na(tempIntraH13), tempIntraH13 := 0]
    } else {
      totalHtable[,tempIntraH13 := 0]
    }
    
    totalHtable[, ':='(tempInterH = tempH - tempIntraH,
                       tempInterH13 = temp13H - tempIntraH13)]
    targettrees <- setkey(targettrees[,k:=1], k)[setkey(weightTable[,k:=1], k),
                                                      nomatch = NA, allow.cartesian = TRUE]
    targettrees[,':='(k = NULL, distanceWeight = NULL, sizeWeight = NULL)]
    totalHtable <- setkey(totalHtable, PlotNumber, competitionName)[setkey(targettrees[,.(PlotNumber, competitionName, TreeNumber, 
                                                                         toDistance)],
                                                          PlotNumber, competitionName), nomatch = 0]
    totalHtable[, overLapArea := 2*(maxRadius^2)*acos(0.5*toDistance/maxRadius)-0.5*toDistance*((4*(maxRadius^2)-toDistance^2)^0.5)]
    HTable <- totalHtable[, .(PlotNumber, TreeNumber, Year = yeardata$Year[1], competitionName, 
                                    H = tempH*500/overLapArea,
                                    IntraH = tempIntraH*500/overLapArea,
                                    InterH = tempInterH*500/overLapArea)]
    totalHtable[,':='(areaRatio = (500/3)/overLapArea)]
    verifyTable <- totalHtable[,.(PlotNumber, TreeNumber, Year = yeardata$Year[1],
                                  competitionName, areaRatio, correctedH = temp13H/areaRatio,
                                  obsH = tempH, correctedIntraH = tempIntraH13/areaRatio,
                                  obsIntraH = tempIntraH, correctedInterH = tempInterH13/areaRatio,
                                  obsInterH = tempInterH)]
                                     
    rm(totalHtable)
    if(i == 1){
      output <- HTable
      verifyTableoutput <- verifyTable
    } else {
      output <- rbind(output, HTable)
      verifyTableoutput <- rbind(verifyTableoutput, verifyTable)
    }
  }
  output <- unique(output, by = c("PlotNumber", "TreeNumber", "Year", "competitionName"))
  competitionNames <- unlist(lapply(weightTable$competitionName, function(x) paste(c("H_", "IntraH_", "InterH_"), x, sep = "")))
  output <- dcast(output, PlotNumber+TreeNumber+Year~competitionName, 
                   value.var = c("H", "IntraH", "InterH"))
  setcolorder(output, c("PlotNumber", "TreeNumber", "Year", competitionNames))
  verifyTableoutput <- unique(verifyTableoutput, by = c("PlotNumber", "TreeNumber", "Year", "competitionName"))
  verifyTableoutput <- dcast(verifyTableoutput, PlotNumber+TreeNumber+Year+areaRatio~competitionName, 
                  value.var = c("correctedH", "obsH", "correctedIntraH", "obsIntraH",
                                "correctedInterH", "obsInterH"))
  competitionNames <- unlist(lapply(weightTable$competitionName, function(x) paste(c("correctedH_", "obsH_", "correctedIntraH_", 
                                                                                     "obsIntraH_", "correctedInterH_", "obsInterH_"), x, sep = "")))
  setcolorder(verifyTableoutput, c("PlotNumber", "TreeNumber", "Year", "areaRatio", competitionNames))
  if(!(testing)){
    return(list(output = output))
  } else {
    return(list(output = output, verifyOutput = verifyTableoutput))
  }
}

