#' @importFrom terra ncell values xyFromCell cellFromXY
getInitialPosition <- function(nbAnimals,initialPosition,spatialCovs,beta){

  spatialcovnames <- names(spatialCovs)

  if(missing(initialPosition)){
    message("   Randomly drawing initial positions from UD...")
    UD <- getUD(spatialCovs, beta=beta,log=TRUE)
    initPos <- matrix(sample(terra::ncell(UD[[1]]),nbAnimals,replace=FALSE,prob=exp(terra::values(UD[[1]]))/sum(exp(terra::values(UD[[1]])))),
                      1,nbAnimals,byrow=TRUE)
    initialPosition <- t(mapply(function(x) terra::xyFromCell(UD,initPos[,x]),1:nbAnimals,SIMPLIFY = FALSE))
  } else {
    if(is.list(initialPosition)){
      if(length(initialPosition)!=nbAnimals) stop("initialPosition must be a list of length ",nbAnimals)
      for(i in 1:nbAnimals){
        if(length(initialPosition[[i]])!=2 | !is.numeric(initialPosition[[i]]) | any(!is.finite(initialPosition[[i]]))) stop("each element of initialPosition must be a finite numeric vector of length 2")
      }
    } else {
      if(length(initialPosition)!=2 | !is.numeric(initialPosition) | any(!is.finite(initialPosition))) stop("initialPosition must be a finite numeric vector of length 2")
      tmpPos<-initialPosition
      initialPosition<-vector('list',nbAnimals)
      for(i in 1:nbAnimals){
        initialPosition[[i]]<-tmpPos
      }
    }
    for(i in 1:length(initialPosition)){
      for(j in 1:length(spatialCovs)){
        if(is.na(terra::cellFromXY(spatialCovs[[j]],matrix(initialPosition[[i]],ncol=2)))) stop("initialPosition for individual ",i," is not within the spatial extent of the ",spatialcovnames[j]," raster")
      }
    }
  }
  initialPosition <- do.call(rbind,initialPosition)
  return(initialPosition)
}
