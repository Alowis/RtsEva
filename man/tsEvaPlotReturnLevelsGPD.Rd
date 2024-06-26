% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{tsEvaPlotReturnLevelsGPD}
\alias{tsEvaPlotReturnLevelsGPD}
\title{tsEvaPlotReturnLevelsGPD}
\usage{
tsEvaPlotReturnLevelsGPD(
  epsilon,
  sigma,
  threshold,
  epsilonStdErr,
  sigmaStdErr,
  thresholdStdErr,
  nPeaks,
  timeHorizonInYears,
  rlvmax,
  tstamps,
  trans,
  ...
)
}
\arguments{
\item{epsilon}{The shape parameter of the GPD.}

\item{sigma}{The scale parameter of the GPD.}

\item{threshold}{The threshold parameter of the GPD.}

\item{epsilonStdErr}{The standard error of the shape parameter.}

\item{sigmaStdErr}{The standard error of the scale parameter.}

\item{thresholdStdErr}{The standard error of the threshold parameter.}

\item{nPeaks}{The number of peaks used in the GPD estimation.}

\item{timeHorizonInYears}{The time horizon in years for the GPD estimation.}

\item{rlvmax}{A data frame containing the return levels of annual maxima.}

\item{tstamps}{The title for the plot.}

\item{trans}{The transformation type for the return levels.}

\item{...}{Additional arguments to be passed to the function.}
}
\value{
A ggplot object representing the plot of return levels.
}
\description{
\code{tsEvaPlotReturnLevelsGPD} is a function that plots the return levels
using the Generalized Pareto Distribution (GPD).
}
\examples{
# Define the required function arguments
epsilon <- 0.2
sigma <- 0.5
threshold <- 10
epsilonStdErr <- 0.05
sigmaStdErr <- 0.05
thresholdStdErr <- 0.1
rlvmax <- data.frame(
  haz.RP = c(2, 5, 10, 20, 50, 100, 200, 500, 1000),
  Idt = as.POSIXct(as.Date("2000-01-01") + round(runif(9, 0, 21 * 365.25)),
    origin = "1970-01-01"
  ),
  QNS = c(10, 12, 13, 13.2, 14, 15.7, 16, 16.2, 18)
)
tstamps <- "Example Timestamps"
trans <- "ori"
nPeaks=70
SampleTimeHorizon=70
# Call the function with the defined arguments
result <- tsEvaPlotReturnLevelsGPD(
  epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr,nPeaks,
  SampleTimeHorizon,rlvmax, tstamps, trans
)
# Plot the result
result
}
\seealso{
\code{\link{tsEvaComputeReturnLevelsGPD}}
\code{\link{tsEvaPlotReturnLevelsGPDFromAnalysisObj}}
}
