#' @title Construct species distribution maps models
#' 
#' @description Constructs mean, standard deviation, and quantile (0.025, 0.5 and 0.975) maps for models calculated using \code{\link{uniSpace}} and \code{\link{ppSpace}} 
#' 
#' @param modelSpace An object of class \code{ppSpace} or of class \code{uniSpace}.
#' @param dims A vector of length 2 defining the number of pixels to use as rows and columns to define the map.
#' @param sPoly A spatial polygon to isolate the region of interest. If none is given, a map is drawn for the entire region covered by the mesh. 
#' @param sample Logical. Whether to sample from the posterior distribution of each cell using \code{INLA}'s \code{inla.posterior.sample}. Returns a raster stack of each sample along with other types specified. Default to \code{FALSE}. Currently ignored.
#' @param nsamples Integer. Number of samples to draw from the posterior. Ignored if \code{sample = FALSE}.
#' 
#' @importFrom INLA inla.mesh.projector
#' @importFrom INLA inla.mesh.project
#' @importFrom INLA inla.stack.index
#' @importFrom INLA inla.posterior.sample
#' @importFrom raster raster
#' @importFrom raster mask
#' @importFrom raster xmin
#' @importFrom raster ymin
#' @importFrom raster xmax
#' @importFrom raster ymax
#' @importFrom terra rast
#' @importFrom terra ext
#'
#' @keywords hplot
#' 
#' @export
#' 
mapSpace <- function(modelSpace, dims, sPoly = NULL, sample = FALSE, nsamples = 100){
  
  ### Names to extract and/or assign to layers
  valsPred<-c("mean", "sd", "0.025quant", "0.5quant", "0.975quant","mode")
  valsLink<-c("mean", "sd", "0.025quant", "0.5quant", "0.975quant","mode")
  valsSpat<-c("mean","sd")
  if(sample){
    valsSamp<-paste0("sample",formatC(1:nsamples,width=nchar(nsamples),flag=0))
  }else{
    valsSamp<-NULL
  }
  
  ### Define map basis
  if(is.null(sPoly)){
    mapBasis <- inla.mesh.projector(attributes(modelSpace)$mesh,
                                    dims = dims,
                                    crs = attributes(modelSpace)$mesh$crs$crs)
  }else{
    mapBasis <- inla.mesh.projector(attributes(modelSpace)$mesh,
                                    dims = dims,
                                    xlim = c(xmin(sPoly), xmax(sPoly)),
                                    ylim = c(ymin(sPoly), ymax(sPoly)),
                                    crs = attributes(modelSpace)$mesh$crs$crs)
  }
  
  

  ### Project spatial field values  
  mapSpat <- as.matrix(inla.mesh.project(mapBasis, 
                                           modelSpace$summary.random[['i']][,valsSpat]))

  
  ### Find the mesh edges on which predictions should be made
  ID <- inla.stack.index(attributes(modelSpace)$Stack, tag="pred")$data
    
  ### Project predicted values
  mapPred <- as.matrix(inla.mesh.project(mapBasis, 
                                           modelSpace$summary.fitted.values[ID,valsPred]))
  ### Project predicted values
  mapLink <- as.matrix(inla.mesh.project(mapBasis, 
                                         modelSpace$summary.linear.predictor[ID,valsLink]))
  
  ### Sample from posterior and project sampled values  
  if(sample){
    class(modelSpace)<-"inla"
    samps<-inla.posterior.sample(nsamples,modelSpace)
    samps<-lapply(samps,function(i){
      i$latent
    })
    samps<-do.call("cbind",samps)
    vals<-samps[grep("i:",row.names(samps)),]
    vals<-samps[1:nrow(vals),] # the second set from nrow(vals)+1 to 2*nrow(vals) is the same
    mapSamp <- as.matrix(inla.mesh.project(mapBasis,vals))
    mat<-cbind(mapPred,mapLink,mapSpat,mapSamp)
  }else{
    mat<-cbind(mapPred,mapLink,mapSpat)
  }
    
    
  ### Rearrange and reorder values to put them in a raster stack
  # Probably not the most efficient way to do that...
  a<-array(rev(as.vector(mat)),dim=c(dims[1],dims[2],ncol(mat)))
  a<-apply(a,3,t,simplify=FALSE)
  a<-lapply(a,function(i){
    i[,ncol(i):1]
  })
  a<-simplify2array(a)
  mapRaster<-rast(a[,,dim(a)[3]:1]) # uses terra for now
  ext(mapRaster)<-c(xmin = min(mapBasis$x), xmax = max(mapBasis$x),ymin = min(mapBasis$y), ymax = max(mapBasis$y))
  names(mapRaster)<-c(valsPred,paste0("link",valsLink),paste0("space",valsSpat),valsSamp)
  mapRaster
}