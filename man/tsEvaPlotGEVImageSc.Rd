% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{tsEvaPlotGEVImageSc}
\alias{tsEvaPlotGEVImageSc}
\title{tsEvaPlotGEVImageSc}
\usage{
tsEvaPlotGEVImageSc(
  Y,
  timeStamps,
  serix,
  epsilon,
  sigma,
  mu,
  returnPeriodInDts,
  maxObs,
  trans,
  varargin
)
}
\arguments{
\item{Y}{A vector of extreme values.}

\item{timeStamps}{A vector of timestamps corresponding to the extreme values.}

\item{serix}{The y-value at which to draw a horizontal line on the plot.}

\item{epsilon}{A numeric value representing the shape parameter of the GEV distribution.}

\item{sigma}{A vector of scale parameters corresponding to the timestamps.}

\item{mu}{A vector of location parameters corresponding to the timestamps.}

\item{returnPeriodInDts}{The return period in decimal time steps.}

\item{maxObs}{A data frame containing the maximum observations.}

\item{trans}{A character string indicating the transformation for the plot.
Possible values are "rev" (reverse), inv" (inverse),
lninv (log of inverse) and "ori"(original).}

\item{varargin}{Additional arguments to customize the plot.}
}
\value{
A ggplot object representing the GEV plot with a raster image.
}
\description{
\code{tsEvaPlotGEVImageSc}is a function that generates a plot of the
Generalized Extreme Value (GEV) distribution with evolving parameters
using the provided data.
}
\seealso{
\code{\link{tsEvaPlotGEVImageScFromAnalysisObj}}
}
