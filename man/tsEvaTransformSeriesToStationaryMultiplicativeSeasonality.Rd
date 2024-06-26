% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trend_functions.R
\name{tsEvaTransformSeriesToStationaryMultiplicativeSeasonality}
\alias{tsEvaTransformSeriesToStationaryMultiplicativeSeasonality}
\title{tsEvaTransformSeriesToStationaryMultiplicativeSeasonality}
\usage{
tsEvaTransformSeriesToStationaryMultiplicativeSeasonality(
  timeStamps,
  series,
  timeWindow,
  seasonalityVar = TRUE
)
}
\arguments{
\item{timeStamps}{A vector of timestamps for the time series data.}

\item{series}{A vector of the time series data.}

\item{timeWindow}{The size of the moving window used for trend estimation.}

\item{seasonalityVar}{A logical value indicating whether to consider
a time varying seasonality (30 years moving average)  or a static
seasonal cycle in the transformation. Default is TRUE.}
}
\value{
A list containing the transformed data and various statistics and errors.
\describe{
\item{\code{runningStatsMulteplicity}}{The size of the moving window used for trend estimation}
\item{\code{stationarySeries}}{The transformed stationary series}
\item{\code{trendSeries}}{The trend component of the transformed series}
\item{\code{trendSeriesNonSeasonal}}{The trend component of the original series without seasonality}
\item{\code{stdDevSeries}}{The standard deviation component of the transformed series}
\item{\code{stdDevSeriesNonSeasonal}}{The standard deviation component of the original series without seasonality}
\item{\code{trendNonSeasonalError}}{The error on the non-seasonal trend component}
\item{\code{stdDevNonSeasonalError}}{The error on the non-seasonal standard deviation component}
\item{\code{trendSeasonalError}}{The error on the seasonal trend component}
\item{\code{stdDevSeasonalError}}{The error on the seasonal standard deviation component}
\item{\code{trendError}}{The overall error on the trend component}
\item{\code{stdDevError}}{The overall error on the standard deviation component}
\item{\code{Regime}}{The estimated regime of the trend seasonality}
\item{\code{timeStamps}}{The input timestamps}
\item{\code{nonStatSeries}}{The original non-stationary series}
\item{\code{statSer3Mom}}{The third moment of the transformed stationary series}
\item{\code{statSer4Mom}}{The fourth moment of the transformed stationary series}
}
}
\description{
This function decomposes a time series into a season-dependent trend and a
season-dependent standard deviation.It performs a transformation from
non-stationary to stationary.
}
\section{transformation non stationary -> stationary}{
transformation stationary -> non stationary
y(t) = stdDev(t)*ssn_stdDev(t)*x(t) + trend(t) + ssn_trend(t)
trasfData.trendSeries = trend(t) + ssn_trend(t)
trasfData.stdDevSeries = stdDev(t)*ssn_stdDev(t)
}

\examples{
timeAndSeries <- ArdecheStMartin
timeStamps <- ArdecheStMartin[,1]
series <- ArdecheStMartin[,2]
#select only the 5 latest years
yrs <- as.integer(format(timeStamps, "\%Y"))
tokeep <- which(yrs>=2015)
timeStamps <- timeStamps[tokeep]
series <- series[tokeep]
timeWindow <- 365 # 1 year
TrendTh <- NA
result <- tsEvaTransformSeriesToStationaryMultiplicativeSeasonality(timeStamps,
series, timeWindow,seasonalityVar=FALSE)
plot(result$trendSeries)
}
