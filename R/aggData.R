#' @importFrom deldir deldir
#' @importFrom deldir tile.list
#' @importFrom sp Polygon
#' @importFrom sp Polygons
#' @importFrom sp SpatialPolygons
#' @importFrom sp over
#' 
aggData <- function(xyt, meshSpace, meshTime=NULL, meshDual){
  
  #================
  ### Basic objects
  #================
  nSpaceEdges <- meshSpace$n
  
  if(!is.null(meshTime)){
    nTimeEdges <- meshTime$n
  }
  
  #============================================
  ### Aggregate data points over space and time
  #============================================
  
  o <- st_intersects(meshDual, xyt)
  spaceAgg <- data.frame(space = 1:nrow(meshDual), Freq = lengths(o))

  
  ### Organise temporal data for aggregation
  if(is.null(meshTime)){
    res <- spaceAgg
  }else{
    timeEdges <- sort(c(meshTime$loc[c(1,nTimeEdges)],
                        meshTime$loc[2:nTimeEdges-1]/2 + 
                          meshTime$loc[2:nTimeEdges]/2))
    
    timeAgg <- factor(findInterval(xyt$t, timeEdges), 
                      levels = 1:nTimeEdges)
    
    ### Aggregate over space and time
    spaceTimeAgg <- as.data.frame(table(rep(factor(spaceAgg$space,levels=1:nrow(meshDual)), spaceAgg$Freq), timeAgg))
    spaceTimeAgg <- apply(spaceTimeAgg, 2, 
                          function(x) as.integer(as.character(x)))  
  
    ### Return
    res <- as.data.frame(spaceTimeAgg)
    colnames(res)[1:2] <- c("space", "time")
  }
  
  return(res)
}