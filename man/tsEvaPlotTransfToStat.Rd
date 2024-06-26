% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{tsEvaPlotTransfToStat}
\alias{tsEvaPlotTransfToStat}
\title{tsEvaPlotTransfToStat}
\usage{
tsEvaPlotTransfToStat(
  timeStamps,
  statSeries,
  srsmean,
  stdDev,
  st3mom,
  st4mom,
  varargin
)
}
\arguments{
\item{timeStamps}{A vector of time stamps for the data points.}

\item{statSeries}{A vector of the main time series data.}

\item{srsmean}{A vector of the mean values for each time stamp.}

\item{stdDev}{A vector of the standard deviation values for each time stamp.}

\item{st3mom}{A vector of the third moment values for each time stamp.}

\item{st4mom}{A vector of the fourth moment values for each time stamp.}

\item{varargin}{Additional optional arguments to customize the plot.}
}
\value{
A ggplot object representing the line plot.
}
\description{
\code{tsEvaPlotTransfToStat}is a function that creates a
line plot of time series data along with statistical measures.
}
\seealso{
\code{\link{tsEvaPlotTransfToStatFromAnalysisObj}}
}
