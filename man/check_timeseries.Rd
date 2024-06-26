% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_timeseries}
\alias{check_timeseries}
\title{Check if all years in a time series are present}
\usage{
check_timeseries(timeseries, yro)
}
\arguments{
\item{timeseries}{A time series object.}

\item{yro}{A vector specifying the start and end years.}
}
\value{
A logical value indicating whether all years in the time series are present.
}
\description{
This function checks if all years specified in a given time series are present.
}
\examples{
ts_data <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2004-12-31"), by = "year")
check_timeseries(ts_data, c(2000, 2004))
# Output: TRUE

check_timeseries(ts_data, c(2000, 2005))
# Output: FALSE

}
