#' @title Univariate spatial model
#' @name uniSpace
#'
#' @description Spatial model for presence-absence data using INLA. This function is essentially a specialized wrapper over \code{inla}
#'
#' @param formula A formula that relates the response and some (or all) of the explanatory variables \code{X}. The name of the response variable is the name of the variable in \code{sPoints}.
#' @param sPoints A \code{SpatialPointsDataFrame} with the response variable.
#' @param explanaMesh An object of class \code{explanaMesh}
#' @param offset A character string defining the explanatory variable in \code{explanaMesh$X} to use as offset.
#' @param family A character string describing the error distribution to be used when constructing the model. 
#' @param link A character string describing the link function to be used when constructing the model.
#' @param smooth A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.
#' @param prior.range A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).
#' @param prior.sigma A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).
#' @param \dots Arguments passed to \code{inla}
#'
#' @details 
#'
#' The underlying model used by this function is a generalized spatial linear model using the cloglog link function. The idea to use the cloglog link function instead of another link function (e.g. logit or probit) is that it is more flexible.
#'
#' @return 
#' 
#' An object of class \code{uniSpace} that includes a model output, which is the model output of INLA.
#' 
#' In addition, it includes a series of attributes:
#' 
#'	  \item{\code{formula}}{The formula used to construct the model}
#'	  \item{\code{sPoints}}{A \code{SpatialPointDataFrame} object that includes the sample location and associated data of the modelled species.}
#'	  \item{\code{XEst}}{A matrix with all the explanatory variables used to construct the model. If there were factors in the original set of explanatory variables \code{X}, in \code{XEst}, they were decomposed into dummy variables. The values in \code{XEst} are the one from the sampled location.}
#'	  \item{\code{XPred}}{A matrix with all the explanatory variables used to construct the model. If there were factors in the original set of explanatory variables \code{X}, in \code{XPred}, they were decomposed into dummy variables. The values in \code{XPred} were gathered at the mesh edges.}
#'	  \item{\code{mesh}}{An object of class \code{inla.mesh}. It is the mesh used to construct the model.}
#'	  \item{\code{Stack}}{An object of class \code{inla.data.stack}. It is a stack object constructed internally.}
#'
#' @importFrom sp coordinates
#' @importFrom raster rasterFromXYZ
#' @importFrom sp SpatialPoints
#' @importFrom raster extract
#' @importFrom INLA inla.spde2.pcmatern
#' @importFrom INLA inla.spde.make.A
#' @importFrom INLA inla.stack
#' @importFrom INLA inla.stack.data
#' @importFrom INLA inla.stack.A
#' @importFrom INLA inla
#' @importFrom stats model.matrix
#' @importFrom stats model.frame
#'
#' @export
#'
#' @keywords models
uniSpace <- function(formula, 
                     sPoints,  
                     explanaMesh,
                     offset = NULL, 
                     family = NULL,
                     link = NULL,
                     smooth = 2,
                     prior.range = c(0.05, 0.01),
                     prior.sigma = c(1, 0.01),
                     ...){
  #============
  # Basic check
  #============
  if(is.null(family) | is.null(link)){
    stop("Either 'family' or 'link' (or both) need to be specified")
  }
  
  if(!is.null(offset)){
    explanForm <- as.character(formula)[3]
    if(length(grep(offset,explanForm)) > 0){
      stop("'offset' is also in 'formula'")
    }
  }
  
  #==============
  # Basic objects
  #==============
  nsmpl <- length(sPoints)
  nameresp <- all.vars(formula[[2]])
  nEdges <- explanaMesh$meshSpace$n

  #============
  # Define SPDE
  #============
  SPDE <- inla.spde2.pcmatern(mesh=explanaMesh$meshSpace,
                              alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.sigma)

  #========================
  # Define response objects
  #========================
  # Coordinates use for the estimation
  xy <- coordinates(sPoints)

  #=========================
  # Define covariate objects
  #=========================
  # Organize data into a data.frame
  refData <- as.data.frame(values(explanaMesh$X))
  
  Xfactor <- unlist(lapply(explanaMesh$Xmesh, is.factor))
  if(any(Xfactor)){
    for(i in which(Xfactor)){
      refData[,i] <- as.factor(refData[,i])
      levels(refData[,i]) <- levels(explanaMesh$Xmesh[,i])
    }
  }
  
  refData <- data.frame(y = 1, refData)
  colnames(refData)[1] <- nameresp

  # Organize X so that it follows the formula
  Xorg <- model.matrix(formula,model.frame(formula,
                                           data = refData,
                                           na.action = NULL))[,-1,drop=FALSE]

  # Construct a brick out of Xorg
  xyXorg <- cbind(coordinates(explanaMesh$X),Xorg)
  Xbrick <- rasterFromXYZ(xyXorg)
  if(nlayers(Xbrick)==1){
    Xbrick <- brick(Xbrick)  
  }
  names(Xbrick) <- colnames(Xorg)
  
  if(!is.null(offset)){
    Xbrick <- addLayer(Xbrick, explanaMesh$X[[offset]])
  }
  
  # Extract covariate values for model estimation
  locEst <- SpatialPoints(coords = xy)
  XEst <- extract(Xbrick, locEst)
  XEst <- as.data.frame(cbind(Intercept = 1, XEst))
  
  # Extract covariate values for model prediction
  locPred <- SpatialPoints(coords = explanaMesh$meshSpace$loc[,1:2])

  prefData <- data.frame(y = 1, explanaMesh$Xmesh)
  colnames(prefData)[1] <- nameresp
  
  XPred <- model.matrix(formula,model.frame(formula,
                                            data = prefData,
                                            na.action = NULL))[,-1,drop=FALSE]
  
  if(!is.null(offset)){
    XPred <- cbind(XPred, explanaMesh$Xmesh[,offset])
    colnames(XPred)[ncol(XPred)] <- offset
  }
  
  XPred <- as.data.frame(cbind(Intercept = 1, XPred))
  
  #===================
  # Construct A matrix
  #===================
  # For estimation
  AEst <- inla.spde.make.A(mesh = explanaMesh$meshSpace, loc = xy)

  # For prediction
  APred<-inla.spde.make.A(mesh = explanaMesh$meshSpace)

  #==================================
  # Build stack object for estimation
  #==================================
  # For estimation
  resp <- unlist(sPoints[,nameresp]@data)
  names(resp) <- NULL
  resp <- list(y = resp)
  names(resp) <- nameresp

  StackEst <- inla.stack(data = resp, A = list(AEst, 1),
                         effects = list(list(i = 1:nEdges),
                                        XEst),
                         tag = "est")

  # For prediction
  respNA <- list(y = NA)
  names(respNA) <- nameresp

  StackPred <- inla.stack(data = respNA, A = list(APred, 1),
                          effects = list(list(i = 1:nEdges),
                                         XPred),
                          tag = "pred")

  # Combine both stack objects
  Stack <- inla.stack(StackEst, StackPred)

  #============
  # Build model
  #============
  if(is.null(offset)){
    forEffect <- paste(colnames(XPred),collapse = "+")
    
    formule <- formula(paste(nameresp ,
                             "~ 0 +", forEffect,
                             "+ f(i, model=SPDE)"))
  }else{
    forEffect <- paste(colnames(XPred)[colnames(XPred) != offset],
                       collapse = "+")
    
    formule <- formula(paste(nameresp ,
                             "~ 0 +", forEffect,
                             "+offset(",offset,") + f(i, model=SPDE)"))
  }

  model <- inla(formule, family = family,
                data = inla.stack.data(Stack),
                control.family = list(link = link),
                control.predictor = list(A = inla.stack.A(Stack),
                                         link = 1), ...)
  
  nameRes <- names(model)
  #=============
  # Return model
  #=============
  # add a series of attributes to the result object
  attributes(model) <- list(formula = formula,
                            sPoints = sPoints,
                            XEst = XEst,
                            XPred = XPred,
                            meshSpace = explanaMesh$meshSpace,
                            Stack = Stack)
  names(model) <- nameRes
  
  class(model) <- "uniSpace"

  return(model)
}

