#' @title Prepare data for model construction
#' @name explanaMesh
#' 
#' @description Prepare and organise the data for model prediction. This function is meant to be used as a starting point for all spatial and spatiotemporal model in this package. This function also gathers explanatory variables values associated to the edges within and outside the area of interest.
#' 
#' @param sPoly An sf polygon defining the sampling region.
#' @param meshSpace An \code{inla.mesh} object for the spatial part of the model.
#' @param meshTime An \code{inla.mesh} object for the temporal part of the model. Default is NULL
#' @param X A \code{\link{terra}} raster that includes the explanatory variables to consider for the analysis. Default is NULL.
#' @param verbose Logical. Whether or not to print on the screen (20 times) the number of edges that have been assigned values.
#' 
#' @details 
#' 
#' The time it takes to run this function is directly related to the number of vertex (points) in the mesh; more specifically the ones outside the region of interest.
#' 
#' Also, this function checks the crs of \code{sPoly}, \code{meshSpace} and \code{X} to make sure they match. Note that all objects are required to have a crs.
#' 
#' For spatiotemporal models, \code{meshTime} needs to be defined. However, it is not a necessary argument because \code{explanaMesh} is also meant to be used to build spatial model without any temporal components.
#' 
#' When building spatiotemporal models, it is possible that explanatory variables change through time. For this reason, \code{X} should be a \code{\link[terra]{SpatRaster}} where layers associated to a particular environment variable that changes through time have the same name but are differentiated by sequential numbers separated by a point. For example, they could have the following names : "env.1", "env.2", "env.3", etc. Time associated to the environment is defined as a z-value in the raster, using \code{\link[terra]{time}}. For all analyses to work properly, it is essential that time be formatted in a \code{\link{POSIXct}} class, otherwise the function will send an error.
#' 
#' If there are no explanatory variable available, it is still possible to estimate the hyper-parameters of a spatial or a spatio-temporal model through the function in the mapSpecies R package.
#' 
#' @return 
#' 
#' An object of class \code{explanaMesh} that includes a list of all the objects used as arguments (except for the verbose call) and \code{Xmesh} the values of the explanatory variables for all edges of the mesh.
#' 
#' @importFrom sf st_as_sf
#' @importFrom sf st_crs
#' @importFrom terra extract
#' @importFrom terra xyFromCell
#' @importFrom terra ncell
#' @importFrom terra cellFromXY
#' @importFrom terra vect
#' @importFrom terra time
#' @importFrom terra ext
#' @importFrom terra xmin
#' @importFrom terra xmax
#' @importFrom terra ymin
#' @importFrom terra ymax
#' @importFrom stats dist
#'
#' @keywords manip
#'
#' @export
#' 
explanaMesh <- function(sPoly, meshSpace, meshTime = NULL, X = NULL, verbose = TRUE){
  #=================
  # Check projection
  #=================
  
  # Extract projection of objects
  projSp <- st_crs(sPoly)
  if(!is.null(X)){
    projX <- st_crs(X)
  }
  projMesh <- st_crs(meshSpace$crs)
  
  if(is.null(X)){
    proj <- list(projSp, projMesh)
  }else{
    proj <- list(projSp, projMesh, projX)
  }
  
  proja<-sapply(proj,function(i){
    i == proj[[1]]
  })
  
  if(!all(proja)){
    stop("sPoly, meshSpace and X do not all have the same crs")
  }
  
  #===================
  # Temporal component
  #===================
  if(!is.null(meshTime)){
    # Check if meshTime is an inla.mesh.1d object
    if(class(meshTime) != "inla.mesh.1d"){
      stop("'meshTime' needs to be an object of class 'inla.mesh.1d'")
    }
    
    # Check for time value
    if(!is.null(X)){
      zval <- time(X)
      if(!any(class(zval) == "POSIXct")){
        stop("'X' needs to have time values with a POSIXct class.")
      }
      
      # Check layers names
      splitXName <- strsplit(names(X), "\\.")
      nSplitXName <- unlist(lapply(splitXName, length))
      
      if(!all(nSplitXName == 2)){
        stop("Some of the variables do not have the structure 'variable.number'")
      }
    }
  }
  
  #============================================================
  # If the edges are not within the raster boundaries associate
  # them with a close by pixel
  #
  # Some of this code was borrows from : 
  # https://stackoverflow.com/questions/26652629/extracting-a-value-from-a-raster-for-a-specific-point-based-on-the-closest-cell
  #============================================================
  if(!is.null(X)){

    # Extract the points
    pts <- meshSpace$loc[,1:2]
    
    # Find the quadrant around the raster in which the pixel is located
    quad <- colSums(sapply(pts[, 1], '>', ext(X)[1:2])) * 3 + 
      colSums(sapply(pts[, 2], '>', ext(X)[3:4]))
    
    Quad <- split(as.data.frame(pts), quad)
    
    # Find new coordinates (by quadrants)
    
    new.pts <- lapply(names(Quad), function(x) {
      switch(x, 
             '0' = matrix(xyFromCell(X, ncell(X) - ncol(X) + 1)[rep(1, nrow(Quad[[x]])), ], ncol = 2),
             '1' = matrix(xyFromCell(X, cellFromXY(X, cbind(xmin(X), Quad[[x]][, 2]))), ncol = 2),
             '2' = matrix(xyFromCell(X, 1)[rep(1, nrow(Quad[[x]])), ], ncol = 2),
             '3' = matrix(xyFromCell(X, cellFromXY(X, cbind(Quad[[x]][, 1], ymin(X)))), ncol = 2),
             '4' = matrix(unlist(Quad[[x]]), ncol = 2),
             '5' = matrix(xyFromCell(X, cellFromXY(X, cbind(Quad[[x]][, 1], ymax(X)))), ncol = 2),
             '6' = matrix(xyFromCell(X, ncell(X))[rep(1, nrow(Quad[[x]])), ], ncol = 2),
             '7' = matrix(xyFromCell(X, cellFromXY(X, cbind(xmax(X), Quad[[x]][, 2]))), ncol = 2),
             '8' = matrix(xyFromCell(X, ncol(X))[rep(1, nrow(Quad[[x]])), ], ncol = 2)
      )
    }
    )
    
    # Reorganize to a data.frame 
    new.pts <- unsplit(mapply(function(x, y) {
      row.names(x) <- row.names(y)
      as.data.frame(x)
    }, 
    new.pts, Quad, SIMPLIFY=FALSE), quad)
    
    # Construct SpatialPoints from meshSpace edges
    loc <- st_as_sf(new.pts, coords=1:2,
                    crs = st_crs(sPoly))
    
    # Extract values from X at loc
    locVal <- extract(X, vect(loc), ID = FALSE) |> as.matrix()
    
    #==================================================
    # If there are some elements in locVal that are NAs
    #==================================================
    # Find locVal with NAs
    locValNA <- sort(unique(which(is.na(locVal), arr.ind = TRUE)[,1]))
    nlocValNA <- length(locValNA)
    
    if(nlocValNA > 0){
      loc <- loc[locValNA, ]
      for(i in 1:nlocValNA){
        extX <- ext(X)
        extDist <- dist(matrix(c(extX$xmin,extX$xmax,
                                 extX$ymin,extX$ymax),nrow = 2))
        buffer <- extDist/10
        locExtract <- extract(X,
                              vect(loc[i,]),
                              buffer = buffer,
                              fun = mean,
                              ID = FALSE)
        locExtract <- as.matrix(locExtract)
        while(any(is.na(locExtract))){
          locExtract <- extract(X, vect(loc[i,]), 
                                buffer = buffer,
                                fun = mean,
                                ID = FALSE)
          locExtract <- as.matrix(locExtract)
          buffer <- buffer + buffer
        }
        
        if(!is.null(locExtract)){
          locVal[locValNA[i],] <- locExtract
        }
        
        locExtract <- NULL
        
        if(verbose){
          if(any(i == round(seq(0, nlocValNA, length = 20)[-c(1,20)]))){ # 20 is the number of times things are printed
            print(paste(i," out of ",nlocValNA," edges with NAs considered",
                        sep = ""))
          }
        }
      }
    }
    
    # Convert to data.frame
    locVal <- as.data.frame(locVal)
    
    # Check for factors in X and convert locVal accordingly
    Xfactor <- is.factor(X)
    
    # Not sure what the following does. Computes the mean level from integer representation and than finds the mean level?
    if(any(Xfactor)){
      for(i in which(Xfactor)){
        locVal[,i] <- round(locVal[,i])
        locVal[,i] <- as.factor(locVal[,i])
        levels(locVal[,i]) <- levels(X)[[i]][,2]
      }
    }
    
    # return
    if(is.null(meshTime)){
      results <- list(sPoly = sPoly, meshSpace = meshSpace, meshTime = NULL, X = X, Xmesh = locVal)
    }else{
      results <- list(sPoly = sPoly, meshSpace = meshSpace, meshTime = meshTime, X = X, Xmesh = locVal)
    }
  }else{
    # return
    if(is.null(meshTime)){
      results <- list(sPoly = sPoly, meshSpace = meshSpace, meshTime = NULL, X = NULL, Xmesh = NULL)
    }else{
      results <- list(sPoly = sPoly, meshSpace = meshSpace, meshTime = meshTime, X = NULL, Xmesh = NULL)
    }
  }
  
  class(results) <- "explanaMesh"
  
  return(results)
}

