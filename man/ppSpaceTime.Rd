% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppSpaceTime.R
\name{ppSpaceTime}
\alias{ppSpaceTime}
\title{Spatio-temporal point process model}
\usage{
ppSpaceTime(
  formula,
  ST,
  ppWeight,
  explanaMesh,
  meshTime,
  timeRes,
  offset = NULL,
  smooth = 2,
  prior.range = c(0.05, 0.01),
  prior.sigma = c(1, 0.01),
  prior.pccor = c(0.7, 0.7),
  many = FALSE,
  ...
)
}
\arguments{
\item{formula}{A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.}

\item{ST}{An object of class \code{\link[spacetime]{STF}}*, \code{\link[spacetime]{STI}}* or \code{\link[spacetime]{STS}}*.}

\item{ppWeight}{An object of class \code{\link{ppWeight}}}

\item{explanaMesh}{An object of class \code{dataPrep}}

\item{meshTime}{An object of class \code{\link{inla.mesh.1d}}}

\item{timeRes}{A character string defining the temporal resolution to use to perform the analysis. The character string choices are given in the Details section of \code{\link{DateTimeClasses}} help.}

\item{offset}{A character string defining the explanatory variable in \code{explanaMesh$X} to use as offset.}

\item{smooth}{A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.}

\item{prior.range}{A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).}

\item{prior.sigma}{A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).}

\item{prior.pccor}{A vector of length 2, with (cor, Pcor) specifying that P(cor > cor_0) = p_cor, where cor is the temporal autocorrelation. Default is c(0.7, 0.7).}

\item{many}{Logical. Whether the data in \code{sPoints} is large or not. See details. Default is \code{FALSE}.}

\item{\dots}{Arguments passed to \code{inla}}
}
\value{
A list including an \code{\link{inla}} object, an\code{\link{inla.stack}} object, an \code{\link{inla.mesh.2d}} object for the spatial component of the model and an \code{\link{inla.mesh.1d}} object for the temporal component of the model.
}
\description{
Spatio-temporal point process model using INLA. This function is essentially a sophisticated wrapper over \code{inla}
}
\details{
If the argument \code{many = TRUE}, the estimation and the prediction will be carried out solely at the mesh edges, whereas when \code{many = FALSE} the estimation will be carried out at the mesh edges and at the sampled location. When the number of samples is very large (e.g. tens of thousands of samples or more) using \code{many = TRUE} can be much more computationally efficient. However, there is a precision trade-off. When \code{many = TRUE}, each sample is associated to an edge and the model is constructed using the number of samples associated to an edge as an importance value. In doing so, some precision is lost at the expense of speed.
}
\keyword{models}
