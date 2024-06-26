% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{tsEvaPlotReturnLevelsGEV}
\alias{tsEvaPlotReturnLevelsGEV}
\title{tsEvaPlotReturnLevelsGEV}
\usage{
tsEvaPlotReturnLevelsGEV(
  epsilon,
  sigma,
  mu,
  epsilonStdErr,
  sigmaStdErr,
  muStdErr,
  rlvmax,
  tstamps,
  trans,
  ...
)
}
\arguments{
\item{epsilon}{The shape parameter of the GEV distribution.}

\item{sigma}{The scale parameter of the GEV distribution.}

\item{mu}{The location parameter of the GEV distribution.}

\item{epsilonStdErr}{The standard error of the shape parameter.}

\item{sigmaStdErr}{The standard error of the scale parameter.}

\item{muStdErr}{The standard error of the location parameter.}

\item{rlvmax}{A data frame containing the return levels of annual maxima.}

\item{tstamps}{The title for the plot.}

\item{trans}{The transformation used to fit the EVD, either "ori" (original)
or "rev" (reverse). "inv" and "lninv" are also available
but in development phase.}

\item{...}{Additional arguments to be passed to the function.}
}
\value{
A ggplot object representing the plot of return levels.
}
\description{
\code{tsEvaPlotReturnLevelsGEV} is a function that plots the return levels
using the Generalized Extreme Value (GEV) distribution.
}
\examples{
# Define the required function arguments
epsilon <- 0.2
sigma <- 0.5
mu <- 10
epsilonStdErr <- 0.05
sigmaStdErr <- 0.05
muStdErr <- 0.1
rlvmax <- data.frame(
  haz.RP = c(2, 5, 10, 20, 50, 100, 200, 500, 1000),
  Idt = as.POSIXct(as.Date("2000-01-01") + round(runif(9, 0, 21 * 365.25)),
    origin = "1970-01-01"
  ),
  QNS = c(10, 12, 13, 13.2, 14, 15.7, 16, 16.2, 18)
)
tstamps <- "Example Timestamps"
trans <- "ori"
# Call the function with the defined arguments
result <- tsEvaPlotReturnLevelsGEV(
  epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr,
  rlvmax, tstamps, trans
)

# Plot the result
result
}
\seealso{
\code{\link{tsEvaComputeReturnLevelsGEV}}
\code{\link{tsEvaPlotReturnLevelsGEVFromAnalysisObj}}
}
