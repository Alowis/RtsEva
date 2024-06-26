% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eva_functions.R
\name{tsEVstatistics}
\alias{tsEVstatistics}
\title{tsEVstatistics}
\usage{
tsEVstatistics(
  pointData,
  alphaCI = 0.95,
  gevMaxima = "annual",
  gevType = "GEV",
  evdType = c("GEV", "GPD"),
  shape_bnd = c(-0.5, 1)
)
}
\arguments{
\item{pointData}{A list containing the dataset of extreme values. It should include the following components:
\describe{
\item{\code{annualMax}}{A vector of annual maximum values}
\item{\code{annualMaxDate}}{A vector of dates corresponding to the annual maximum values}
\item{\code{monthlyMax}}{A matrix of monthly maximum values}
}}

\item{alphaCI}{The confidence level for the confidence intervals of the parameter estimates. Default is 0.95.}

\item{gevMaxima}{The type of maxima to use for GEV fitting. Can be either 'annual' or 'monthly'. Default is 'annual'.}

\item{gevType}{The type of GEV distribution to use. Can be either 'GEV', 'Gumbel'. Default is 'GEV'.}

\item{evdType}{The types of extreme value distributions to calculate. Can be a combination of 'GEV' and 'GPD'. Default is c('GEV', 'GPD').}

\item{shape_bnd}{The lower and upper bounds for the shape parameter of the GEV distribution. Default is c(-0.5, 1).}
}
\value{
A list containing the following components:
\describe{
\item{\code{EVmeta}}{A list containing metadata about the analysis.
It includes Tr, A vector of return periods for which
return levels are calculated}
\item{\code{EVdata}}{A list containing the calculated statistics and return levels. It includes the following components:}
\describe{
\item{\code{GEVstat}}{A list containing the GEV statistics and return levels:}
\describe{
\item{\code{method}}{The method used for fitting the GEV distribution.}
\item{\code{values}}{A vector of return levels calculated using the GEV distribution.}
\item{\code{parameters}}{A vector of parameter estimates for the GEV distribution.}
\item{\code{paramCIs}}{A matrix of confidence intervals for the parameter estimates.}
}
}
\item{\code{GPDstat}}{list containing the GPD statistics and return levels:}
\describe{
\item{\code{method}}{The method used for fitting the GPD distribution}
\item{\code{values}}{A vector of return levels calculated using the GPD distribution}
\item{\code{parameters}}{A vector of parameter estimates for the GPD distribution}
\item{\code{paramCIs}}{A matrix of confidence intervals for the parameter estimates}
}
\item{\code{isValid}}{A logical value indicating whether the analysis was performed or not.}
}
}
\description{
\code{tsEvstatistics} is a function that calculates the Generalized Extreme
Value (GEV) and Generalized Pareto Distribution (GPD) statistics
and return levels for a given dataset of extreme values.
}
\examples{
# Create a sample dataset
data <- ArdecheStMartin
colnames(data) <- c("Date", "Value")
yrs <- as.integer(format(data$Date, "\%Y"))
tokeep <- which(yrs>=2015)
data <- data[tokeep,]
pointData <- tsEvaSampleData(data, meanEventsPerYear=3, minEventsPerYear=0,
minPeakDistanceInDays=7, "high")
result <- tsEVstatistics(pointData)
result$EVdata$GEVstat$values
result$EVdata$GPDstat$values
}
