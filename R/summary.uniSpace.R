#' @title Model summary for object of class \code{uniSpace}
#' @method summary uniSpace
#' 
#' @description Extracts a summary of the model parameters estimated with the \code{\link{uniSpace}} function.
#' 
#' @param object An object of class \code{\link{uniSpace}}
#' @param \dots Additional arguments affecting the summary produced.
#' 
#' @details
#' 
#' This summary function extract a series of estimated model parameters and statistics (such as the DIC or the WAIC if they have been calculated) and organises them so that they can be readily interpreted. All the information is available in the \code{\link{uniSpace}} object but it is either difficult to find or badly formatted. This function simplifies the access to this information.
#' 
#' If some non-standard model parameters or statistics have been computed, e.g. through the \code{\link[INLA]{control.compute}} argument, they have to be access directly from the \code{\link{uniSpace}} object.
#' 
#' @return 
#' 
#' A list with the following structure 
#' 
#'    \item{\code{coefficients}}{All the estimated parameters associated to the explanatory variables. Associated to the estimated coefficients, the standard deviation, the 0.025, 0.5 and 0.975 quantile, the mode and the Kullback-Leibler divergence (\code{kld}). The \code{kld} describes the difference between the standard Gaussian approximation and the INLA approximation to the marginal posterior densities. Smaller \code{kld} are better.}
#'    \item{\code{spaceCoefficients}}{Range and standard deviation estimation, standard deviation as well as the 0.025, 0.5 and 0.975 quantile and the mode for the spatial component of the model. Note that \code{prior.range} and \code{prior.sigma} in \code{\link{uniSpace}} are the priors associated to these two estimated hyper parameters.}
#'    \item{\code{nEffectiveParam}}{Number of effective parameters in the model. Because this is a hierarchical model with a spatial component and that information is shared across parameters, the values given are the effective (or independent) number of parameters used in the model. The estimated value is an average and as such also includes a standard deviation.  }
#'    \item{\code{nEquivalentRep}}{Number of equivalent replicates. This value is obtained from dividing the total number of samples by the number of effective parameters in the model. }
#'    \item{\code{runningTime}}{Time to perform each step of the INLA analysis. There are some additional preparation step in the \code{\link{uniSpace}} function, but these should not take more than a few seconds to run at most and should thus account only for a minimal fraction of the time related to computation.}
#'    \item{\code{marginalLike}}{Gaussian and INLA log marginal-likelihood. This argument may be \code{NULL} if it was defined as such by the user through the \code{\link[INLA]{control.compute}} INLA argument.}
#'    \item{\code{WAIC}}{The full model Watanabe-Akaike Information Criterion. This argument may be \code{NULL} if it was defined as such by the user through the \code{\link[INLA]{control.compute}} INLA argument. To access the more details WAIC information calculated through INLA, they are available in the \code{\link{uniSpace}} object.}
#'    \item{\code{DIC}}{The full model Deviance Information Criterion. This argument may be \code{NULL} if it was defined as such by the user through the \code{\link[INLA]{control.compute}} INLA argument. To access the more details DIC information calculated through INLA, they are available in the \code{\link{uniSpace}} object.}
#'
#' @keywords manip
#' 
#' @export
summary.uniSpace <- function(object,...){
  # Extract the parameters associated to the explanatory variables
  fixed <- object$summary.fixed
  
  # Extract the parameters associated to the spatial components
  space <- object$summary.hyperpar
  if(!is.null(space)){
    rownames(space) <- gsub(" for i","",rownames(space)) 
  }
  
  # Processing time
  time <- object$cpu.used[4]
  
  # Result object
  res <- list(coefficients = fixed,
              spaceCoefficients = space,
              nEffectiveParam = object$neffp[1:2,],
              nEquivalentRep = object$neffp[3,],
              runningTime = time,
              marginalLike = object$mlik,
              WAIC = object$waic$waic,
              DIC = object$dic$dic)

  class(res) <- "summary.uniSpace"
  
  return(res)
}
