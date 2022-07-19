#' @title Calculate weight for the integration of the point process
#' 
#' @description This function calculates the weight that should be given to the 
#' 
#' @param sPoly A \code{SpatialPolygons*} object.
#' @param mesh An \code{inla.mesh} object
#' 
#' @return 
#' 
#' An object of class \code{ppWeight} which returns a vector of weights. 
#' 
#' @importFrom sf st_as_sf
#' @importFrom sf st_intersects
#' @importFrom sf st_intersection
#' @importFrom sf st_area
#'
#' @export
#' 
#' @keywords models
#' 
ppWeight <- function(sPoly, mesh){
  #-------------------------------------------------
  ### Construct a dual mesh from the triangular mesh
  #-------------------------------------------------
  dmesh <- inla.mesh.dual(mesh)

  #--------------------------------------------------------------
  ### Find the intersection between the polygons in the dual mesh
  ### and the location domain
  #--------------------------------------------------------------

  ### Calculate weight
  weight <- numeric(nrow(dmesh))
  overlaps <- st_intersects(dmesh, st_as_sf(sPoly))
  cuts <- st_intersection(dmesh, st_as_sf(sPoly))
  w <- which(as.logical(sapply(overlaps, length)))
  areas <- as.numeric(st_area(cuts))
  if(length(w)!=length(areas)){
    stop("Number of resulting polygons different from the number of touching cells")
  }
  weight[w] <- areas
  
  ### Check to make sure there are integration points with 0 weight
  if(all(weight > 0)){
    stop("There needs to be some weights that are of 0")
  }
  
  ### Return mesh
  attributes(weight) <- list(mesh=mesh,dmesh=as(dmesh,"Spatial"))
  
  class(weight) <- "ppWeight"
  
  return(weight)
}