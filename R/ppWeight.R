#' @title Calculate weight for the integration of the point process
#' 
#' @description This function calculates the weight that should be given to the 
#' 
#' @param sPoly An sf polygon object.
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
  
  # An intersect and a within are used to first find polygons that will 
  # be cut to reduce the duration of the intersection (about twice as fast)
  
  ### Calculate weight
  dmesh$id <- 1:nrow(dmesh)
  weight <- numeric(nrow(dmesh))
  overlaps <- lengths(st_intersects(dmesh, sPoly))
  within <- lengths(st_within(dmesh, sPoly))
  o <- overlaps > 0L & !within > 0L
  suppressWarnings(
    cuts <- st_intersection(sPoly, dmesh[o, ])
  )
  dmeshcuts <- rbind(cuts, dmesh[within > 0L, ])
  dmeshcuts <- dmeshcuts[order(dmeshcuts$id), ]
  w <- which(as.logical(overlaps))
  areas <- as.numeric(st_area(dmeshcuts))
  if(length(w)!=length(areas)){
    stop("Number of resulting polygons different from the number of touching cells")
  }
  weight[w] <- areas
  
  ### Check to make sure there are integration points with 0 weight
  if(all(weight > 0)){
    stop("There needs to be some weights that are of 0")
  }
  
  ### Return mesh
  attributes(weight) <- list(mesh = mesh, dmesh = dmesh, dmeshcuts = dmeshcuts)
  
  class(weight) <- "ppWeight"
  
  return(weight)
}