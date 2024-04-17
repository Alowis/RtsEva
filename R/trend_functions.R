# Detrending main functions-------------------

#' Transform Time Series to Stationary Trend Only
#'
#' This function is the original detrending function implemented in Mentaschi et al.(2016).
#' It takes a time series and transforms it into a stationary one.
#' It computes the trend as a running average of the time series,
#' the slowly varying amplitude as its standard deviation, and other statistical measures.
#'
#' @param timeStamps A vector of time stamps for the time series.
#' @param series The original time series.
#' @param timeWindow The size of the time window used for detrending.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{runningStatsMulteplicity}: The running statistics multiplicity.
#'   \item \code{stationarySeries}: The transformed stationary series.
#'   \item \code{trendSeries}: The trend series.
#'   \item \code{trendSeriesNonSeasonal}: The non-seasonal trend series.
#'   \item \code{trendError}: The error on the trend.
#'   \item \code{stdDevSeries}: The slowly varying standard deviation series.
#'   \item \code{stdDevSeriesNonSeasonal}: The non-seasonal slowly varying standard deviation series.
#'   \item \code{stdDevError}: The error on the standard deviation.
#'   \item \code{timeStamps}: The time stamps.
#'   \item \code{nonStatSeries}: The original non-stationary series.
#'   \item \code{statSer3Mom}: The third moment of the transformed stationary series.
#'   \item \code{statSer4Mom}: The fourth moment of the transformed stationary series.
#' }
#'
#' @examples
#' timeStamps <- seq(as.Date("2022-01-01"), as.Date("2022-01-10"), by = "day")
#' series <- rnorm(length(timeStamps))
#' timeWindow <- 3
#' result <- tsEvaTransformSeriesToStationaryTrendOnly(timeStamps, series, timeWindow)
#'
#' @export
tsEvaTransformSeriesToStationaryTrendOnly <- function(timeStamps, series, timeWindow) {
  cat("computing the trend ...\n")
  rs <- tsEvaDetrendTimeSeries(timeStamps, series, timeWindow)
  nRunMn <- rs@nRunMn
  cat("computing the slowly varying standard deviation ...\n")
  varianceSeries <- tsEvaNanRunningVariance(rs@detrendSeries, nRunMn)
  # further smoothing
  varianceSeries <- tsEvaNanRunningMean(varianceSeries, ceiling(nRunMn / 2))

  stdDevSeries <- varianceSeries^.5
  statSeries <- rs@detrendSeries / stdDevSeries
  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  kurstosis <- moments::kurtosis(statSeries)
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))

  # N is the size of each sample used to compute the average
  N <- nRunMn
  # the error on the trend is computed as the error on the average:
  #  error(average) = stdDev/sqrt(N)
  trendError <- mean(stdDevSeries) / N^.5

  # variance(stdDev) ~ 2 stdDev^4 / (n - 1)
  # stdDevError is approximated as constant
  avgStdDev <- mean(stdDevSeries)
  S <- 2

  # computation of the error on the standard deviation explained in
  # Mentaschi et al 2016
  stdDevError <- avgStdDev * (2 * S^2 / N^3)^(1 / 4)

  trasfData <- list()
  trasfData$runningStatsMulteplicity <- rs@nRunMn
  trasfData$stationarySeries <- statSeries
  trasfData$trendSeries <- rs@trendSeries
  trasfData$trendSeriesNonSeasonal <- rs@trendSeries
  trasfData$trendError <- trendError
  trasfData$stdDevSeries <- stdDevSeries
  trasfData$stdDevSeriesNonSeasonal <- stdDevSeries
  trasfData$stdDevError <- stdDevError * rep(1, length(stdDevSeries))
  trasfData$timeStamps <- timeStamps
  trasfData$nonStatSeries <- series
  trasfData$statSer3Mom <- statSer3Mom
  trasfData$statSer4Mom <- statSer4Mom
  return(trasfData)
}

#' Transform Time Series to Stationary using percentiles to compute the slowly varying amplitude of the distribution.
#'
#' This function transforms a time series to a stationary ones using a moving average as the trend and a running percentiles
#' to represent the slowly varying amplitude of the distribution
#'
#' @param timeStamps A vector of time stamps for the time series.
#' @param series The original time series.
#' @param timeWindow The size of the moving window used for detrending.
#' @param percentile The percentile value used to compute the extreme trend of the stationary series.
#' @param ... Additional arguments.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{runningStatsMulteplicity}{The running statistics multiplicity.}
#'   \item{stationarySeries}{The transformed stationary trend only series.}
#'   \item{trendSeries}{The trend series.}
#'   \item{trendSeriesNonSeasonal}{The non-seasonal trend series.}
#'   \item{trendError}{The error on the trend.}
#'   \item{stdDevSeries}{The standard deviation series.}
#'   \item{stdDevSeriesNonSeasonal}{The non-seasonal standard deviation series.}
#'   \item{stdDevError}{The error on the standard deviation.}
#'   \item{timeStamps}{The time stamps.}
#'   \item{nonStatSeries}{The original non-stationary series.}
#'   \item{statSer3Mom}{The running mean of the third moment of the stationary series.}
#'   \item{statSer4Mom}{The running mean of the fourth moment of the stationary series.}
#' }
#'
#' @examples
#' timeStamps <- seq(as.Date("2022-01-01"), as.Date("2022-01-31"), by = "day")
#' series <- rnorm(length(timeStamps))
#' result <- tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile(timeStamps, series, 7, 95, 0.5)
#'
#' @import stats
#' @export
tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile <- function(timeStamps, series, timeWindow, percentile, ...) {

  # Removing trend from the series
  cat("\ncomputing the trend on extremes...\n")
  qtes <- quantile(series, percentile / 100, na.rm = T)
  rs <- tsEvaDetrendTimeSeries(timeStamps, series, timeWindow)
  meantrend <- mean(series)
  trendseriesf <- rs@trendSeries
  detrendSeries <- rs@detrendSeries
  nRunMn <- rs@nRunMn

  # Compute the extreme trend of the stationary series
  qtes <- quantile(detrendSeries, percentile / 100, na.rm = T)
  detrendseriep <- detrendSeries
  detrendseriep[which(detrendseriep < qtes)] <- NA
  percentileSeries <- tsEvaNanRunningMean(detrendseriep, nRunMn)
  percentileSeries[which(is.na(percentileSeries))] <- min(percentileSeries, na.rm = T)

  # further smoothing
  percentileSeries <- tsEvaNanRunningMean(percentileSeries, ceiling(nRunMn / 2))
  if (min(percentileSeries, na.rm = T) < 0) percentileSeries <- percentileSeries - min(percentileSeries)

  meanperc <- mean(detrendseriep, na.rm = T)
  stdDev <- sd(detrendSeries, na.rm = T)

  # real std
  varianceSeries <- tsEvaNanRunningVariance(detrendSeries, nRunMn)
  # further smoothing
  varianceSeries <- tsEvaNanRunningMean(varianceSeries, ceiling(nRunMn / 2))
  stdDevSeries1 <- varianceSeries^.5

  stdDevSeries <- percentileSeries / meanperc * stdDev


  avgStdDev <- mean(stdDevSeries)
  S <- 2
  # N is the size of each sample used to compute the average
  N <- nRunMn
  # computation of the error on the standard deviation explained in
  # Mentaschi et al 2016
  stdDevError <- avgStdDev * (2 * S^2 / N^3)^(1 / 4)

  # Compute running statistics of the stationary series
  statSeries <- detrendSeries / stdDevSeries
  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))

  # the error on the trend is computed as the error on the average:
  trendError <- mean(stdDevSeries) / N^.5


  # Store the results in a list
  trasfData <- list(
    runningStatsMulteplicity = nRunMn,
    stationarySeries = statSeries,
    trendSeries = trendseriesf,
    trendSeriesNonSeasonal = NULL,
    trendError = trendError,
    stdDevSeries = stdDevSeries,
    stdDevSeriesNonSeasonal = NULL,
    stdDevError = stdDevError * rep(1, length(stdDevSeries)),
    timeStamps = timeStamps,
    nonStatSeries = series,
    statSer3Mom = statSer3Mom,
    statSer4Mom = statSer4Mom
  )
}


#' Transform Time Series to Stationary only using extremes
#'
#' This function transforms a time series to a stationary one by focusing on extremes.
#' The trend and slowly varying amplitude are computed on values above a
#' threshold defined by the user.
#'
#' @param timeStamps A vector of time stamps corresponding to the observations in the series.
#' @param series A vector of the time series data.
#' @param timeWindow The size of the time window used for detrending.
#' @param TrendTh The threshold for fitting the trend on the means above a certain quantile. Default is 0.5.
#' @param ... Additional arguments to be passed to other functions.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{runningStatsMulteplicity}: The multiplicity of running statistics.
#'     \item \code{stationarySeries}: The stationary series after removing the trend.
#'     \item \code{trendSeries}: The trend component of the series.
#'     \item \code{trendSeriesNonSeasonal}: NULL (not used).
#'     \item \code{trendError}: The error on the trend component.
#'     \item \code{stdDevSeries}: The standard deviation series.
#'     \item \code{stdDevSeriesNonSeasonal}: NULL (not used).
#'     \item \code{stdDevError}: The error on the standard deviation series.
#'     \item \code{timeStamps}: The time stamps.
#'     \item \code{nonStatSeries}: The original non-stationary series.
#'     \item \code{statSer3Mom}: The running mean of the third moment of the stationary series.
#'     \item \code{statSer4Mom}: The running mean of the fourth moment of the stationary series.
#'   }
#'
#' @import stats
#' @examples
#' timeStamps <- seq(as.Date("2022-01-01"), as.Date("2022-01-31"), by = "day")
#' series <- rnorm(length(timeStamps))
#' transformedData <- tsEvaTransformSeriesToStationaryPeakTrend(timeStamps, series, timeWindow = 7)
#'
#' @seealso [tsEvaFindTrendThreshold()]
#' @export
tsEvaTransformSeriesToStationaryPeakTrend <- function(timeStamps, series, timeWindow, TrendTh, ...) {
  # Removing trend from the series
  cat("\ncomputing the trend on extremes...\n")

  # fit a trend on the means above a threshold
  # The threshold can be specified as an input
  # If not, default is 0.5
  if (is.na(TrendTh)) {
    TrendTh <- 0.5
  }
  print(paste0("trend threshold= ", TrendTh))
  qd <- quantile(series, TrendTh, na.rm = T)
  serieb <- series
  serieb[which(serieb < qd)] <- NA


  rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow)

  # adding the threshold to avoid values on the zero
  detrendSeries <- series - rs@trendSeries
  detrendSerie1 <- serieb - rs@trendSeries
  qd2 <- min(detrendSerie1, na.rm = T)
  nRunMn <- rs@nRunMn

  varianceSeries <- tsEvaNanRunningVariance(detrendSerie1, nRunMn)
  # further smoothing
  varianceSeries <- tsEvaNanRunningMean(varianceSeries, ceiling(nRunMn / 2))

  stdDevSeries1 <- varianceSeries^.5
  stdDevSeries <- stdDevSeries1
  avgStdDev <- mean(stdDevSeries)
  S <- 2

  # N is the size of each sample used to compute the average
  # Counts the real amount of timstep, need to specify original temporal resolution
  N <- timeWindow * 4
  # computation of the error on the standard deviation explained in
  # Mentaschi et al 2016
  stdDevError <- avgStdDev * (2 * S^2 / N^3)^(1 / 4)

  # Compute running statistics of the stationary series
  statSeries <- detrendSeries / stdDevSeries

  xtremS <- statSeries
  xtremS <- na.omit(xtremS)
  allS <- na.omit(serieb)

  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))


  # the error on the trend is computed as the error on the average:
  #  error(average) = stdDev/sqrt(N)
  trendError <- mean(stdDevSeries) / N^.5

  # Store the results in a list
  trasfData <- list(
    runningStatsMulteplicity = nRunMn,
    stationarySeries = statSeries,
    trendSeries = rs@trendSeries,
    trendSeriesNonSeasonal = NULL,
    trendError = trendError,
    stdDevSeries = stdDevSeries,
    stdDevSeriesNonSeasonal = NULL,
    stdDevError = stdDevError * rep(1, length(stdDevSeries)),
    timeStamps = timeStamps,
    nonStatSeries = series,
    statSer3Mom = statSer3Mom,
    statSer4Mom = statSer4Mom
  )

  return(trasfData)
}


#' Transform Series to Stationary with Multiplicative Seasonality
#'
#' This function decomposes a time series into a season-dependent trend and a season-dependent standard deviation.
#' It performs a transformation from non-stationary to stationary.
#'
#' @param timeStamps A vector of timestamps for the time series data.
#' @param series A vector of the time series data.
#' @param timeWindow The size of the moving window used for trend estimation.
#' @param seasonalityVar A logical value indicating whether to consider seasonality in the transformation. Default is TRUE.
#' @param ... Additional arguments to be passed to other functions.
#'
#' @return A list containing the transformed data and various statistics and errors.
#' \describe{
#'   \item{runningStatsMulteplicity}{The size of the moving window used for trend estimation.}
#'   \item{stationarySeries}{The transformed stationary series.}
#'   \item{trendSeries}{The trend component of the transformed series.}
#'   \item{trendSeriesNonSeasonal}{The trend component of the original series without seasonality.}
#'   \item{stdDevSeries}{The standard deviation component of the transformed series.}
#'   \item{stdDevSeriesNonSeasonal}{The standard deviation component of the original series without seasonality.}
#'   \item{trendNonSeasonalError}{The error on the non-seasonal trend component.}
#'   \item{stdDevNonSeasonalError}{The error on the non-seasonal standard deviation component.}
#'   \item{trendSeasonalError}{The error on the seasonal trend component.}
#'   \item{stdDevSeasonalError}{The error on the seasonal standard deviation component.}
#'   \item{trendError}{The overall error on the trend component.}
#'   \item{stdDevError}{The overall error on the standard deviation component.}
#'   \item{Regime}{The estimated regime of the trend seasonality.}
#'   \item{timeStamps}{The input timestamps.}
#'   \item{nonStatSeries}{The original non-stationary series.}
#'   \item{statSer3Mom}{The third moment of the transformed stationary series.}
#'   \item{statSer4Mom}{The fourth moment of the transformed stationary series.}
#' }
#'
#' @examples
#' timeStamps <- seq(as.Date("2022-01-01"), as.Date("2022-12-31"), by = "day")
#' series <- sin(2 * pi * seq(0, 1, length.out = length(timeStamps)))
#' transformedData <- tsEvaTransformSeriesToStationaryMultiplicativeSeasonality(timeStamps, series, 30)
#' @export

tsEvaTransformSeriesToStationaryMultiplicativeSeasonality <- function(timeStamps, series, timeWindow, seasonalityVar = T, ...) {
  # this function decomposes the series into a season-dependent trend and a
  # season-dependent standard deviation.
  # The season-dependent standard deviation is given by a seasonal factor
  # multiplied by a slowly varying standard deviation.

  # transformation non stationary -> stationary
  # x(t) = [y(t) - trend(t) - ssn_trend(t)]/[stdDev(t)*ssn_stdDev(t)]
  # transformation stationary -> non stationary
  # y(t) = stdDev(t)*ssn_stdDev(t)*x(t) + trend(t) + ssn_trend(t)

  # trasfData.trendSeries = trend(t) + ssn_trend(t)
  # trasfData.stdDevSeries = stdDev(t)*ssn_stdDev(t)
  svar <- 1
  if (seasonalityVar == T) svar <- 2
  seasonalityTimeWindow <- 2 * 30.4 # 2 months

  cat("computing trend ...\n")
  rs <- tsEvaDetrendTimeSeries(timeStamps, series, timeWindow)
  nRunMn <- rs@nRunMn
  cat("computing trend seasonality ...\n")
  statSeries <- rs@detrendSeries
  trendSeasonality <- tsEstimateAverageSeasonality(timeStamps, statSeries, nRunMn)
  Regime <- trendSeasonality$regime
  trendSeasonality <- trendSeasonality$Seasonality
  statSeries <- statSeries - trendSeasonality[, svar]

  cat("computing slowly varying standard deviation ...\n")
  varianceSeries <- tsEvaNanRunningVariance(statSeries, nRunMn)
  # further smoothing
  varianceSeries <- tsEvaNanRunningMean(varianceSeries, ceiling(nRunMn / 2))

  seasonalVarNRun <- round(nRunMn / timeWindow * seasonalityTimeWindow)
  # seasonalVarSeries is a moving variance computed on a short time
  # window of 1-3 months, able to vary with season.
  cat("computing standard deviation seasonality ...\n")
  seasonalVarSeries <- tsEvaNanRunningVariance(statSeries, seasonalVarNRun)
  seasonalStdDevSeries <- sqrt(seasonalVarSeries / varianceSeries)
  seasonalStdDevSeries <- tsEstimateAverageSeasonality(timeStamps, seasonalStdDevSeries, nRunMn)$Seasonality
  seasonalStdDevSeries <- seasonalStdDevSeries[, svar]
  stdDevSeriesNonSeasonal <- sqrt(varianceSeries)
  statSeries <- statSeries / (stdDevSeriesNonSeasonal * seasonalStdDevSeries)
  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))

  # N is the size of each sample used to compute the average
  N <- nRunMn
  # the error on the trend is computed as the error on the average:
  #  error(average) = stdDev/sqrt(N)
  trendNonSeasonalError <- mean(stdDevSeriesNonSeasonal, na.rm = TRUE) / sqrt(N)

  # computation of the error on the standard deviation explained in
  # Mentaschi et al 2016
  S <- 2
  avgStdDev <- mean(stdDevSeriesNonSeasonal, na.rm = TRUE)
  stdDevNonSeasonalError <- avgStdDev * (2 * S^2 / N^3)^(1 / 4)

  Ntot <- length(series)
  trendSeasonalError <- stdDevNonSeasonalError * sqrt(12 / Ntot + 1 / N)
  stdDevSeasonalError <- seasonalStdDevSeries * (288 / Ntot^2 / N)^(1 / 4)

  trendError <- sqrt(trendNonSeasonalError^2 + trendSeasonalError^2)
  stdDevError <- sqrt((stdDevSeriesNonSeasonal * stdDevSeasonalError)^2 + (seasonalStdDevSeries * stdDevNonSeasonalError)^2)

  trasfData <- list()

  trasfData$runningStatsMulteplicity <- nRunMn
  trasfData$stationarySeries <- statSeries
  trasfData$trendSeries <- rs@trendSeries + trendSeasonality[, svar]
  trasfData$trendSeriesNonSeasonal <- rs@trendSeries
  trasfData$stdDevSeries <- stdDevSeriesNonSeasonal * seasonalStdDevSeries
  trasfData$stdDevSeriesNonSeasonal <- stdDevSeriesNonSeasonal
  trasfData$trendNonSeasonalError <- trendNonSeasonalError
  trasfData$stdDevNonSeasonalError <- stdDevNonSeasonalError
  trasfData$trendSeasonalError <- trendSeasonalError
  trasfData$stdDevSeasonalError <- stdDevSeasonalError
  trasfData$trendError <- trendError
  trasfData$stdDevError <- stdDevError
  trasfData$Regime <- Regime
  trasfData$timeStamps <- timeStamps
  trasfData$nonStatSeries <- series
  trasfData$statSer3Mom <- statSer3Mom
  trasfData$statSer4Mom <- statSer4Mom
  return(trasfData)
}


#' Transform Time Series to Stationary with seasonal trend and amplitude
#' based on extreme percentile
#'
#' This function decomposes a time series into a season-dependent trend and a season-dependent standard deviation.
#' The season-dependent amplitude is given by a seasonal factor multiplied by a slowly varying percentile.
#'
#' @param timeStamps A vector of time stamps for the time series.
#' @param series The original time series.
#' @param timeWindow The length of the moving window used for trend estimation.
#' @param percentile The percentile value used for computing the slowly varying percentile.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{runningStatsMulteplicity}{The size of each sample used to compute the average.}
#'   \item{stationarySeries}{The transformed stationary series.}
#'   \item{trendSeries}{The trend series.}
#'   \item{trendSeriesNonSeasonal}{The non-seasonal trend series.}
#'   \item{stdDevSeries}{The season-dependent standard deviation series.}
#'   \item{stdDevSeriesNonSeasonal}{The non-seasonal standard deviation series.}
#'   \item{trendError}{The error on the trend.}
#'   \item{stdDevError}{The error on the standard deviation.}
#'   \item{statSer3Mom}{The 3rd moment of the transformed stationary series.}
#'   \item{statSer4Mom}{The 4th moment of the transformed stationary series.}
#'   \item{nonStatSeries}{The original non-stationary series.}
#'   \item{Regime}{The regime of the trend seasonality.}
#'   \item{timeStamps}{The time stamps.}
#'   \item{trendNonSeasonalError}{The error on the non-seasonal trend.}
#'   \item{stdDevNonSeasonalError}{The error on the non-seasonal standard deviation.}
#'   \item{trendSeasonalError}{The error on the seasonal trend.}
#'   \item{stdDevSeasonalError}{The error on the seasonal standard deviation.}
#' }
#'
#' @examples
#' timeStamps <- seq(as.Date("2022-01-01"), as.Date("2022-12-31"), by = "day")
#' series <- rnorm(length(timeStamps))
#' timeWindow <- 30
#' percentile <- 90
#' result <- tsEvaTransformSeriesToStatSeasonal_ciPercentile(timeStamps, series, timeWindow, percentile)
#'
#' @references
#' Mentaschi, L., Vousdoukas, M. I., Voukouvalas, E., Dosio, A., & Feyen, L. (2016).
#' Joint effect of global sea level rise and tidal changes on the future
#' flood risk in the Rhine-Meuse delta. Environmental Research Letters, 11(3), 035003.
#'
#' @export
tsEvaTransformSeriesToStatSeasonal_ciPercentile <- function(timeStamps, series, timeWindow, percentile) {
  # this function decomposes the series into a season-dependent trend and a
  # season-dependent standard deviation.
  # The season-dependent standard deviation is given by a seasonal factor
  # multiplied by a slowly varying standard deviation.

  # transformation non stationary -> stationary
  # x(t) = [y(t) - trend(t) - ssn_trend(t)]/[stdDev(t)*ssn_stdDev(t)]
  # transformation stationary -> non stationary
  # y(t) = stdDev(t)*ssn_stdDev(t)*x(t) + trend(t) + ssn_trend(t)

  # trasfData.trendSeries = trend(t) + ssn_trend(t)
  # trasfData.stdDevSeries = stdDev(t)*ssn_stdDev(t)

  seasonalityTimeWindow <- 2 * 30.4 # 2 months

  print("computing trend...")
  rs <- tsEvaDetrendTimeSeries(timeStamps, series, timeWindow, percent = percentile)
  nRunMn <- rs@nRunMn
  print("computing trend seasonality...")
  trendSeasonality <- tsEstimateAverageSeasonality(timeStamps, rs@detrendSeries, nRunMn)
  Regime <- trendSeasonality$regime
  trendSeasonality <- trendSeasonality$Seasonality
  statSeries <- rs@detrendSeries - trendSeasonality[, 1]
  print(paste0("computing the slowly varying ", percentile, "th percentile..."))
  percentileSeries <- tsEvaNanRunningPercentiles(timeStamps, statSeries, nRunMn, percentile)
  #
  seasonalVarNRun <- round(nRunMn / timeWindow * seasonalityTimeWindow)
  # seasonalVarSeries is a moving variance computed on a short time
  # window of 1-3 months, able to vary with season.
  print("computing standard deviation seasonality...")
  seasonalVarSeries <- tsEvaNanRunningVariance(statSeries, seasonalVarNRun)
  seasonalStdDevSeries <- sqrt(seasonalVarSeries / percentileSeries$rnprcnt)
  seasonalStdDevSeries <- tsEstimateAverageSeasonality(timeStamps, seasonalStdDevSeries, nRunMn)
  seasonalStdDevSeries <- seasonalStdDevSeries$Seasonality[, 1]
  stdDevSeriesNonSeasonal <- sqrt(percentileSeries$rnprcnt)

  statSeries <- statSeries / (stdDevSeriesNonSeasonal * seasonalStdDevSeries)
  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))

  # N is the size of each sample used to compute the average
  N <- nRunMn
  # the error on the trend is computed as the error on the average:
  #  error(average) = stdDev/sqrt(N)
  trendNonSeasonalError <- mean(stdDevSeriesNonSeasonal) / N^.5

  # computation of the error on the standard deviation explained in
  # Mentaschi et al 2016
  S <- 2
  avgStdDev <- mean(stdDevSeriesNonSeasonal)
  stdDevNonSeasonalError <- avgStdDev * (2 * S^2 / N^3)^(1 / 4)
  Ntot <- length(series)
  trendSeasonalError <- stdDevNonSeasonalError * sqrt(12 / Ntot + 1 / N)
  stdDevSeasonalError <- seasonalStdDevSeries * (288. / Ntot^2 / N)^(1. / 4.)

  trendError <- sqrt(trendNonSeasonalError^2 + trendSeasonalError^2)
  stdDevError <- sqrt((stdDevSeriesNonSeasonal * stdDevSeasonalError)^2 + (seasonalStdDevSeries * stdDevNonSeasonalError)^2)

  trasfData <- list(
    runningStatsMulteplicity = nRunMn,
    stationarySeries = statSeries,
    trendSeries = rs@trendSeries + trendSeasonality[, 1],
    trendSeriesNonSeasonal = rs@trendSeries,
    stdDevSeries = stdDevSeriesNonSeasonal * seasonalStdDevSeries,
    stdDevSeriesNonSeasonal = stdDevSeriesNonSeasonal,
    trendError = trendError,
    stdDevError = stdDevError,
    statSer3Mom = statSer3Mom,
    statSer4Mom = statSer4Mom,
    nonStatSeries = series,
    Regime = Regime,
    timeStamps = timeStamps
  )

  trasfData$trendNonSeasonalError <- trendNonSeasonalError
  trasfData$stdDevNonSeasonalError <- stdDevNonSeasonalError
  trasfData$trendSeasonalError <- trendSeasonalError
  trasfData$stdDevSeasonalError <- stdDevSeasonalError

  return(trasfData)
}


#' Transform Time Series to Stationary Trend and Change Points
#'
#' This function takes a time series and transforms it into a stationary trend series
#' by removing the trend component and detecting change points. It computes the slowly
#' varying standard deviation and normalizes the stationary series before detecting
#' step changes. It also calculates the error on the trend and standard deviation.
#'
#' @param timeStamps A vector of time stamps corresponding to the data points in the series.
#' @param series The original time series data.
#' @param timeWindow The size of the time window used for detrending.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{runningStatsMulteplicity}: The running statistics multiplicity.
#'   \item \code{stationarySeries}: The transformed stationary series.
#'   \item \code{trendSeries}: The trend series.
#'   \item \code{trendonlySeries}: The trend series without the stationary component.
#'   \item \code{ChpointsSeries2}: The trend component of the change points.
#'   \item \code{changePoints}: The detected change points.
#'   \item \code{trendSeriesNonSeasonal}: The trend series without the seasonal component.
#'   \item \code{trendError}: The error on the trend.
#'   \item \code{stdDevSeries}: The slowly varying standard deviation series.
#'   \item \code{stdDevSeriesNonStep}: The slowly varying standard deviation series without step changes.
#'   \item \code{stdDevError}: The error on the standard deviation.
#'   \item \code{timeStamps}: The time stamps.
#'   \item \code{nonStatSeries}: The original non-stationary series.
#'   \item \code{statSer3Mom}: The running mean of the third moment of the stationary series.
#'   \item \code{statSer4Mom}: The running mean of the fourth moment of the stationary series.
#' }
#'
#' @examples
#' timeStamps <- seq(as.Date("2022-01-01"), as.Date("2022-01-10"), by = "day")
#' series <- rnorm(length(timeStamps))
#' timeWindow <- 5
#' result <- tsEvaTransformSeriesToStationaryTrendAndChangepts(timeStamps, series, timeWindow)
#'
#' @export
tsEvaTransformSeriesToStationaryTrendAndChangepts <- function(timeStamps, series, timeWindow) {
  cat('computing the trend ...\n')
  #ChgPtsRaw=tsEvaChangepts(series,timeWindow/2,timeStamps)
  rs <- tsEvaDetrendTimeSeries(timeStamps, series, timeWindow)
  nRunMn = rs@nRunMn
  cat('computing the slowly varying standard deviation ...\n')
  varianceSeries <- tsEvaNanRunningVariance(rs@detrendSeries, nRunMn)
  #further smoothing
  varianceSeries <- tsEvaNanRunningMean(varianceSeries, ceiling(nRunMn/2))
  #
  stdDevSeries <- varianceSeries^.5
  stdDevSeries_1<-stdDevSeries
  statSeries <- rs@detrendSeries
  statSeries_1=statSeries
  cat('computing the step change detection...\n')
  #I normalize the stationary serie before the step change detection
  #Maybe not the best
  statSeries=statSeries/stdDevSeries
  ChgPts=tsEvaChangepts(statSeries,nRunMn/1.5,timeStamps)

  #static std
  statSD=sd(statSeries)
  #ChgptStdDevSeries <- sqrt(ChgPts$variance)/statSD
  ChgptStdDevSeries=1
  statSeries=statSeries-ChgPts$trend
  # statSeries=statSeries/stdDevSeries
  stdDevSeries=stdDevSeries*ChgptStdDevSeries
  trendSeries=rs@trendSeries+ChgPts$trend*stdDevSeries
  # plot(ChgPts$trend)
  # plot(stdDevSeries_1)
  # plot(trendSeries,type="l")
  # lines(rs@trendSeries,col=2)
  #

  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  kurstosis=moments::kurtosis(statSeries)
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))

  # N is the size of each sample used to compute the average
  N <- nRunMn
  # the error on the trend is computed as the error on the average:
  #  error(average) = stdDev/sqrt(N)
  trendError <- mean(stdDevSeries)/N^.5

  # variance(stdDev) ~ 2 stdDev^4 / (n - 1)

  # stdDevError is approximated as constant
  avgStdDev <- mean(stdDevSeries)
  S <- 2
  # computation of the error on the standard deviation explained in
  # Mentaschi et al 2016
  stdDevError <- avgStdDev*( 2*S^2/N^3 )^(1/4)

  trasfData <- list()
  trasfData$runningStatsMulteplicity <- rs@nRunMn
  trasfData$stationarySeries <- statSeries
  trasfData$trendSeries <- trendSeries
  trasfData$trendonlySeries <- rs@trendSeries
  # trasfData$ChpointsSeries <- ChgPtsRaw$trend
  trasfData$ChpointsSeries2 <- ChgPts$trend
  trasfData$changePoints <- ChgPts$changepoints
  trasfData$trendSeriesNonSeasonal <- trendSeries
  trasfData$trendError <- trendError
  trasfData$stdDevSeries <- stdDevSeries
  trasfData$stdDevSeriesNonStep <- stdDevSeries_1
  trasfData$stdDevError <- stdDevError*rep(1,length(stdDevSeries))
  trasfData$timeStamps <- timeStamps
  trasfData$nonStatSeries <- series
  trasfData$statSer3Mom <- statSer3Mom
  trasfData$statSer4Mom <- statSer4Mom
  return(trasfData)
}

#' Transform Time Series to Stationary Trend and Change Points with Confidence Intervals
#'
#' This function takes a time series and transforms it into a stationary trend series with change points and confidence intervals.
#'
#' @param timeStamps A vector of time stamps corresponding to the observations in the series.
#' @param series The time series data.
#' @param timeWindow The size of the sliding window used for detrending the series.
#' @param percentile The percentile value used for computing the running percentile of the stationary series.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{runningStatsMulteplicity}{The running statistics multiplicity.}
#'   \item{stationarySeries}{The transformed stationary series.}
#'   \item{trendSeries}{The trend series.}
#'   \item{trendonlySeries}{The trend series without the stationary component.}
#'   \item{ChpointsSeries2}{The trend series with change points.}
#'   \item{changePoints}{The detected change points.}
#'   \item{trendSeriesNonSeasonal}{The trend series without the seasonal component.}
#'   \item{trendError}{The error on the trend.}
#'   \item{stdDevSeries}{The standard deviation series.}
#'   \item{stdDevSeriesNonStep}{The standard deviation series without the step change component.}
#'   \item{stdDevError}{The error on the standard deviation.}
#'   \item{timeStamps}{The time stamps.}
#'   \item{nonStatSeries}{The original non-stationary series.}
#'   \item{statSer3Mom}{The running mean of the third moment of the stationary series.}
#'   \item{statSer4Mom}{The running mean of the fourth moment of the stationary series.}
#' }
#'
#' @examples
#' timeStamps <- c(1, 2, 3, 4, 5)
#' series <- c(10, 12, 15, 11, 13)
#' timeWindow <- 2
#' percentile <- 90
#' result <- tsEvaTransformSeriesToStationaryTrendAndChangepts_ciPercentile(timeStamps, series, timeWindow, percentile)
#' print(result)
#'
#' @export
tsEvaTransformSeriesToStationaryTrendAndChangepts_ciPercentile <- function(timeStamps, series, timeWindow, percentile) {
  cat('computing the trend ...\n')
  #ChgPtsRaw=tsEvaChangepts(series,timeWindow/2,timeStamps)
  rs <- tsEvaDetrendTimeSeries(timeStamps, series, timeWindow,percent=percentile)
  nRunMn = rs@nRunMn
  cat('computing the slowly varying standard deviation ...\n')
  # Compute the running percentile of the stationary series
  percentileSeries <- tsEvaNanRunningPercentiles(timeStamps, rs@detrendSeries, nRunMn, percentile)
  if(min(percentileSeries$rnprcnt)<0) percentileSeries$rnprcnt=percentileSeries$rnprcnt-min(percentileSeries$rnprcnt)
  meanperc=mean(percentileSeries$rnprcnt,na.rm=T)

  #here there is a problem
  stdDev = sd(rs@detrendSeries,na.rm=T);
  stdDevSeries <- percentileSeries$rnprcnt/meanperc*stdDev
  stdDevError =percentileSeries$stdError/meanperc*stdDev;
  stdDevSeries_1<-stdDevSeries
  statSeries <- rs@detrendSeries
  statSeries_1=statSeries
  cat('computing the step change detection...\n')
  #I normalize the stationary serie before the step change detection
  #Maybe not the best
  statSeries=statSeries/stdDevSeries
  ChgPts=tsEvaChangepts(statSeries,nRunMn/1.5,timeStamps)

  #static std
  statSD=sd(statSeries)
  #ChgptStdDevSeries <- sqrt(ChgPts$variance)/statSD
  ChgptStdDevSeries=1
  statSeries=statSeries-ChgPts$trend
  # statSeries=statSeries/stdDevSeries
  stdDevSeries=stdDevSeries*ChgptStdDevSeries
  trendSeries=rs@trendSeries+ChgPts$trend*stdDevSeries
  # plot(ChgPts$trend)
  # plot(stdDevSeries_1)
  # plot(trendSeries,type="l")
  # lines(rs@trendSeries,col=2)
  #

  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  kurstosis=moments::kurtosis(statSeries)
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))

  # N is the size of each sample used to compute the average
  N <- nRunMn
  # the error on the trend is computed as the error on the average:
  #  error(average) = stdDev/sqrt(N)
  trendError <- mean(stdDevSeries)/N^.5

  # variance(stdDev) ~ 2 stdDev^4 / (n - 1)

  # stdDevError is approximated as constant
  avgStdDev <- mean(stdDevSeries)
  S <- 2
  # computation of the error on the standard deviation explained in
  # Mentaschi et al 2016
  stdDevError <- avgStdDev*( 2*S^2/N^3 )^(1/4)

  trasfData <- list()
  trasfData$runningStatsMulteplicity <- rs@nRunMn
  trasfData$stationarySeries <- statSeries
  trasfData$trendSeries <- trendSeries
  trasfData$trendonlySeries <- rs@trendSeries
  #trasfData$ChpointsSeries <- ChgPtsRaw$trend
  trasfData$ChpointsSeries2 <- ChgPts$trend
  trasfData$changePoints <- ChgPts$changepoints
  trasfData$trendSeriesNonSeasonal <- trendSeries
  trasfData$trendError <- trendError
  trasfData$stdDevSeries <- stdDevSeries
  trasfData$stdDevSeriesNonStep <- stdDevSeries_1
  trasfData$stdDevError <- stdDevError*rep(1,length(stdDevSeries))
  trasfData$timeStamps <- timeStamps
  trasfData$nonStatSeries <- series
  trasfData$statSer3Mom <- statSer3Mom
  trasfData$statSer4Mom <- statSer4Mom
  return(trasfData)
}

# Helpers for main trend functions---------------

#' Find Trend Threshold
#'
#' This function calculates the optimal trend threshold for a given time series.
#'
#' @param series The time series data.
#' @param timeStamps The timestamps corresponding to the time series data.
#' @param timeWindow The time window for detrending the time series.
#'
#' @import stats
#' @return The trend threshold value.
#'
#' @details This function iterates over different percentiles and calculates the threshold based on each percentile. It then removes data points below the threshold and detrends the time series using the specified time window. The function calculates the correlation between the normalized trend and the time series and stores the correlation coefficient for each percentile. It performs a changepoint analysis to determine if there is a significant change in the correlation coefficients. If a change point is found, the function returns the percentile corresponding to the change point. If no change point is found, the function returns the percentile with the highest correlation coefficient. If there are negative values in the detrended time series, the function returns the percentile with the fewest negative values.
#'
#' @examples
#' series <- c(1, 2, 3, 4, 5)
#' timeStamps <- c("2022-01-01", "2022-01-02", "2022-01-03", "2022-01-04", "2022-01-05")
#' timeWindow <- 3
#' find_trend_threshold(series, timeStamps, timeWindow)
#'
#' @export
tsEvaFindTrendThreshold <- function(series, timeStamps, timeWindow){
  # Initialize variables
  ptn=timeStamps[which(!is.na(timeAndSeries$data))]
  bounds=unique(lubridate::year(ptn))
  nr <- rep(1, length(series))
  sts <- c()
  lnegs=c()
  pctd=c()
  #plot(series)
  #I play with this
  pcts <- seq(0.40, 0.96, by=0.02)
  # Iterate over each percentile
  for (iter in 1:length(pcts)){
    # Calculate the threshold based on the current percentile
    thrsdt <- quantile(series, pcts[iter], na.rm=TRUE)
    # Assign a flag value to missing data
    series_no_na <- series
    series_no_na[which(is.na(series_no_na))] <- -9999
    # Remove data points below the threshold
    serieb <- series_no_na
    timeb=timeStamps
    timeb=timeb[-which(serieb < thrsdt)]
    serieb[which(serieb < thrsdt)] <- NA
    #first I have to test that all years are represented
    checkY=check_timeseries(timeb,bounds)
    if (checkY==FALSE){
      break
    }
    # Detrend the time series
    rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow,fast=T)
    varianceSeries <- tsEvaNanRunningVariance(serieb, rs@nRunMn)
    varianceSeries <- tsEvaNanRunningMean(varianceSeries, ceiling(rs@nRunMn/2))
    norm_trend <- rs@trendSeries/mean(rs@trendSeries, na.rm=TRUE)
    # if(length(which(is.na(varianceSeries)))>0) {
    #   break
    # }
    # if(length(which(is.na(norm_trend)))>0) {
    #   break
    # }
    dtr1=serieb-rs@trendSeries
    lneg=length(which(dtr1<0))
    # Calculate the correlation between the normalized trend and the time series
    stab <- cor(nr, norm_trend, use = 'pairwise.complete.obs')
    # Update nr with the current norm_trend if it's the first iteration
    if(iter == 1) nr <- norm_trend
    # Store the correlation coefficient in sts
    if (lneg>=1) lnegs=c(lnegs,lneg)
    sts <- c(sts, stab)
    pctd=c(pctd,pcts[iter])


  }
  # Set the first correlation coefficient to 1
  sts[1] <- 1
  # Perform a changepoint analysis
  cptm <- try(changepoint::cpt.var(sts, method='AMOC',penalty="BIC", minseglen=3),T)
  #plot(cptm)
  if(inherits(cptm, "try-error")){
    rval=pctd[which.max(sts[-1])]

  }else if(cptm@cpts[1] == length(sts)) {
    print("no change point")
    cptm@cpts[1] <- length(sts)/2
    rval=pctd[cptm@cpts[1]]
  }
  if (sum(lnegs)>1){
    rval=pctd[which.min(lnegs)]
  }
  print(rval)
  return(rval)
}

#' Change point detection in time series
#'
#' This function applies the PELT method for change point detection in a time series.
#' It returns the mean and variance of the series segments between change points.
#'
#' @param series A numeric vector representing the time series.
#' @param timeWindow An integer specifying the minimum length of segments.
#' @param timeStamps A vector of timestamps corresponding to the series data points.
#'
#' @return A list containing:
#' \itemize{
#'   \item{trend}{A numeric vector of the same length as series, with each segment between change points filled with its mean value.}
#'   \item{variance}{A numeric vector of the same length as series, with each segment between change points filled with its variance.}
#'   \item{changepoints}{A vector of timestamps at which change points were detected.}
#' }
#'
#' @import changepoint
#' @examples
#' \dontrun{
#' series <- rnorm(100)
#' timeStamps <- 1:100
#' result <- tsEvaChangepts(series, 10, timeStamps)
#' plot(series, type = "l")
#' points(timeStamps[result$changepoints], result$trend[result$changepoints], col = "red")
#' }
#' @export
tsEvaChangepts <- function(series, timeWindow, timeStamps) {
  if (min(series) < 0) {
    shift <- min(series)
  } else {
    shift <- 0
  }
  series <- series - 1.1 * shift
  cptm <- changepoint::cpt.meanvar(series, penalty = "Manual", method = "PELT", pen.value = 100, minseglen = timeWindow, test.stat = "Gamma", shape = 1)
  cpts <- c(1, cpts(cptm), length(series)) # change point time points\
  cptsmean <- cptm@param.est$scale * cptm@param.est$shape + 1.1 * shift

  meants <- rep(0, length(series))
  varts <- rep(0, length(series))
  for (s in 1:length(cptsmean)) {
    meants[cpts[s]:cpts[s + 1]] <- cptsmean[s]
  }
  changepoint <- timeStamps[cpts[-c(1, length(cpts))]]
  return(list(trend = meants, variance = varts, changepoints = changepoint))
}

#' Initialize Percentiles
#'
#' This function calculates the percentiles and probabilities for a given dataset.
#'
#' @param subsrs The input dataset.
#' @param percentM The percentile for the lower bound.
#' @param percent The percentile for the middle bound.
#' @param percentP The percentile for the upper bound.
#'
#' @return A list containing the calculated percentiles and probabilities.
#' @import stats
#' @seealso [tsEvaNanRunningPercentiles()]
#' @export
#'
#' @examples
#' subsrs <- c(1, 2, 3, 4, 5)
#' initPercentiles(subsrs, 25, 50, 75)
initPercentiles <- function(subsrs, percentM, percent, percentP) {
  probObj <- list()
  probObj$percentM <- percentM
  probObj$percent <- percent
  probObj$percentP <- percentP

  probObj$probM <- percentM / 100
  probObj$prob <- percent / 100
  probObj$probP <- percentP / 100

  probObj$tM <- quantile(subsrs, probs = percentM / 100, na.rm = T)
  probObj$t <- quantile(subsrs, probs = percent / 100, na.rm = T)
  probObj$tP <- quantile(subsrs, probs = percentP / 100, na.rm = T)

  probObj$N <- sum(!is.na(subsrs))
  return(probObj)
}

#' Calculate the running mean of a time series with NaN handling
#'
#' This function calculates the running mean of a time series, taking into account NaN values.
#' It uses a sliding window approach to calculate the mean, where the window size is specified by the user.
#' If the number of non-NaN values within the window is greater than a threshold, the mean is calculated.
#' Otherwise, NaN is returned.
#'
#' @param series The input time series
#' @param windowSize The size of the sliding window
#' @return A vector containing the running mean of the time series
#'
#' @examples
#' series <- c(1, 2, NaN, 4, 5, 6, NaN, 8, 9)
#' windowSize <- 3
#' result <- tsEvaNanRunningMean(series, windowSize)
#' print(result)
#'
#' @export
tsEvaNanRunningMean <- function(series, windowSize) {
  minNThreshold <- 1

  rnmn <- matrix(nrow = length(series), ncol = 1)
  dx <- floor(windowSize / 2)
  l <- length(series)
  sm <- 0
  n <- 0
  smm <- c()
  snne <- c()
  spp <- c()
  for (ii in c(1:l)) {
    minindx <- max(ii - dx, 1)
    maxindx <- min(ii + dx, l)
    if (ii == 1) {
      subsrs <- series[minindx:maxindx]
      sm <- sum(subsrs, na.rm = T)
      n <- sum(!is.na(subsrs))
    } else {
      if (minindx > 1) {
        sprev <- series[minindx - 1]
        if (!is.na(sprev)) {
          sm <- sm - sprev
          n <- n - 1
        }
      }
      if (maxindx < l) {
        snext <- series[maxindx]
        if (!is.na(snext)) {
          sm <- sm + snext
          n <- n + 1
        }
      }
    }
    if (n > minNThreshold) {
      rnmn[ii] <- sm / n
    } else {
      rnmn[ii] <- NaN
    }
  }
  return(as.vector(rnmn))
}

#' Calculate the running variance of a time series with NaN handling
#'
#' This function calculates the running variance of a time series, taking into account NaN values.
#' The series must be zero-averaged before passing it to this function.
#'
#' @param series The time series data.
#' @param windowSize The size of the window used for calculating the running variance.
#'
#' @return A vector containing the running variance values.
#'
#' @examples
#' series <- c(1, 2, 3, NaN, 5, 6, 7)
#' windowSize <- 3
#' tsEvaNanRunningVariance(series, windowSize)
#'
#' @export
tsEvaNanRunningVariance <- function(series, windowSize) {
  # !!! series must be 0 averaged!!

  minNThreshold <- 1

  rnmn <- matrix(nrow = length(series), ncol = 1)
  dx <- ceiling(windowSize / 2)
  l <- length(series)
  smsq <- 0
  n <- 0

  for (ii in c(1:l)) {
    minindx <- max(ii - dx, 1)
    maxindx <- min(ii + dx, l)
    if (ii == 1) {
      subSqSrs <- series[minindx:maxindx]^2
      smsq <- sum(subSqSrs, na.rm = T)
      n <- sum(!is.na(subSqSrs))
    } else {
      if (minindx > 1) {
        sprev <- series[minindx - 1]
        if (!is.na(sprev)) {
          smsq <- max(0, smsq - sprev^2)
          n <- n - 1
        }
      }
      if (maxindx < l) {
        snext <- series[maxindx + 1]
        if (!is.na(snext)) {
          smsq <- smsq + snext^2
          n <- n + 1
        }
      }
    }
    if (n > minNThreshold) {
      rnmn[ii] <- smsq / n
    } else {
      rnmn[ii] <- NaN
    }
    # print(n)
  }
  return(rnmn)
}

#' tsEvaNanRunningStatistics
#'
#' Returns the moving statistical momentums to the forth.
#'
#' @param series The input time series data.
#' @param windowSize The size of the moving window.
#'
#' @return A data frame containing the following running statistics:
#' \itemize{
#'   \item{rnvar}{running variance}
#'   \item{rn3mom}{running third statistical momentum}
#'   \item{rn4mom}{running fourth statistical momentum}
#' }
#'
#' @details This function calculates the running variance, running third statistical momentum, and running fourth statistical momentum for a given time series data using a moving window approach. The window size determines the number of observations used to calculate the statistics at each point.
#'
#' @examples
#' series <- c(1, 2, 3, 4, 5)
#' windowSize <- 3
#' tsEvaNanRunningStatistics(series, windowSize)
#'
#' @export
tsEvaNanRunningStatistics <- function(series, windowSize) {
  # tsEvaNanRunningStatistics: returns the the moving statistical momentums
  # to the forth.
  # rnvar: running variance
  # rn3mom: running third statistical momentum
  # rn4mom: running fourth statistical momentum

  minNThreshold <- 1

  rnmn <- tsEvaNanRunningMean(series, windowSize)
  rnvar <- matrix(nrow = length(series), ncol = 1)
  rn3mom <- matrix(nrow = length(series), ncol = 1)
  rn4mom <- matrix(nrow = length(series), ncol = 1)

  dx <- ceiling(windowSize / 2)
  l <- length(series)
  sm <- 0
  smsq <- 0
  sm3pw <- 0
  sm4pw <- 0
  n <- 0
  for (ii in c(1:l)) {
    # cat(paste0(" ",n))
    minindx <- max(ii - dx, 1)
    maxindx <- min(ii + dx, l)
    if (ii == 1) {
      subsrs <- series[minindx:maxindx]
      subsrsMean <- rnmn[1]
      subSqSrs <- (subsrs - subsrsMean)^2
      tg <- moments::skewness(subsrs)
      kg <- moments::kurtosis(subsrs)
      sub3pwSrs <- (subsrs - subsrsMean)^3
      sub4pwSrs <- (subsrs - subsrsMean)^4
      smsq <- sum(subSqSrs, na.rm = T)
      sm3pw <- sum(sub3pwSrs, na.rm = T)
      sm4pw <- sum(sub4pwSrs, na.rm = T)
      n <- sum(!is.na(subSqSrs))
      tes <- sm3pw / n
    } else {
      if (minindx > 1) {
        sprev <- series[minindx - 1] - rnmn[minindx - 1]
        if (!is.na(sprev)) {
          smsq <- max(0, smsq - sprev^2)
          sm3pw <- sm3pw - sprev^3
          sm4pw <- max(0, sm4pw - sprev^4)
          n <- n - 1
        }
      }
      if (maxindx < l) {
        snext <- series[maxindx + 1] - rnmn[minindx + 1]
        if (!is.na(snext)) {
          sm <- sm + snext
          smsq <- smsq + snext^2
          sm3pw <- sm3pw + snext^3
          sm4pw <- sm4pw + snext^4
          n <- n + 1
        }
      }
    }
    if (n > minNThreshold) {
      rnvar[ii] <- smsq / n
      rn3mom[ii] <- sm3pw / n
      rn4mom[ii] <- sm4pw / n
    } else {
      rnvar[ii] <- NaN
      rn3mom[ii] <- NaN
      rn4mom[ii] <- NaN
    }
  }
  output <- data.frame(rnvar, rn3mom, rn4mom)
  return(output)
}


#' tsEvaNanRunningPercentiles
#'
#' Computes a running percentile for a given series using a window with a specified size.
#'
#' @param timeStamps The timestamps of the series.
#' @param series The input series.
#' @param windowSize The size of the window for the running percentile. Must be greater than or equal to 100.
#' @param percent The percent level to which the percentile is computed.
#' @return A list containing the approximated running percentile (rnprcnt) and the standard error (stdError).
#' @details
#' This function computes a running percentile for a given series using a window with a specified size.
#' The running percentile is computed by interpolating the percentile value for the requested percentage
#' based on the quitting values and incoming values in the window.
#' The function also performs smoothing on the output and calculates the standard error.
#'
#' The function uses the following label parameters:
#' \describe{
#'   \item{percentDelta}{Delta for the computation of a percentile interval around the requested percentage.
#'                      If the windowSize is greater than 2000, percentDelta is set to 1.
#'                      If the windowSize is between 1000 and 2000, percentDelta is set to 2.
#'                      If the windowSize is between 100 and 1000, percentDelta is set to 5.}
#'   \item{nLowLimit}{Minimum number of non-NA elements for a window for percentile computation.}
#' }
#'
#' @examples
#' timeStamps <- c(1, 2, 3, 4, 5)
#' series <- c(10, 20, 30, 40, 50)
#' windowSize <- 3
#' percent <- 90
#' result <- tsEvaNanRunningPercentiles(timeStamps, series, windowSize, percent)
#' print(result$rnprcnt)
#' print(result$stdError)
#'
#' @import stats
#' @export
tsEvaNanRunningPercentiles <- function(timeStamps, series, windowSize, percent) {
  if (windowSize > 2000) {
    percentDelta <- 1
  } else if (windowSize > 1000) {
    percentDelta <- 2
  } else if (windowSize > 100) {
    percentDelta <- 5
  } else {
    stop("window size cannot be less than 100")
  }
  nLowLimit <- 100
  percentP <- percent + percentDelta
  if (percentP > 100) {
    stop(paste0("max percent: ", 100 - percentDelta))
  }
  percentM <- percent - percentDelta
  if (percentM < 0) {
    stop(paste0("min percent: ", percentDelta))
  }

  percentP <- percent + percentDelta
  if (percentP > 100) {
    stop(paste0("max percent: ", 100 - percentDelta))
  }

  rnprcnt0 <- rep(NA, length(series))
  dx <- ceiling(windowSize / 2)
  l <- length(series)
  minindx <- 1
  maxindx <- min(1 + dx, l)
  subsrs <- series[minindx:maxindx]
  probObj <- initPercentiles(subsrs, percentM, percent, percentP)
  rnprcnt0[1] <- probObj$t

  for (ii in 2:l) {
    minindx <- max(ii - dx, 1)
    maxindx <- min(ii + dx, l)
    if (minindx > 1) {
      sprev <- series[minindx - 1]

      # removing element and reviewing probability
      if (!is.na(sprev)) {
        Nold <- probObj$N
        Nnew <- probObj$N - 1

        nle <- as.numeric(sprev < probObj$tM)
        probObj$probM <- (probObj$probM * Nold - nle) / Nnew
        probObj$percentM <- probObj$probM * 100

        nle <- as.numeric(sprev < probObj$t)
        probObj$prob <- (probObj$prob * Nold - nle) / Nnew
        probObj$percent <- probObj$prob * 100

        nle <- as.numeric(sprev < probObj$tP)
        probObj$probP <- (probObj$probP * Nold - nle) / Nnew
        probObj$percentP <- probObj$probP * 100

        probObj$N <- Nnew
      }
    }
    if (maxindx < l) {
      snext <- series[maxindx + 1]

      if (!is.na(snext)) {
        Nold <- probObj$N
        Nnew <- probObj$N + 1

        nle <- as.numeric(snext < probObj$tM)
        probObj$probM <- (probObj$probM * Nold + nle) / Nnew
        probObj$percentM <- probObj$probM * 100

        nle <- as.numeric(snext < probObj$t)
        probObj$prob <- (probObj$prob * Nold + nle) / Nnew
        probObj$percent <- probObj$prob * 100

        nle <- as.numeric(snext < probObj$tP)
        probObj$probP <- (probObj$probP * Nold + nle) / Nnew
        probObj$percentP <- probObj$probP * 100

        probObj$N <- Nnew
      }
    }

    cout1 <- probObj$percentM > percent
    cout2 <- probObj$percentP < percent
    outOfInterval <- cout1 || cout2
    if (outOfInterval) {
      subsrs <- series[minindx:maxindx]
      probObj <- initPercentiles(subsrs, percentM, percent, percentP)
    }
    if (probObj$N > nLowLimit) {
      if (percent == probObj$percentM) {
        prcntii <- probObj$tM
      } else if ((probObj$percentM < percent) && (percent < probObj$percent)) {
        h1 <- probObj$percent - percent
        h2 <- percent - probObj$percentM
        prcntii <- (h1 * probObj$tM + h2 * probObj$t) / (h1 + h2)
      } else if (percent == probObj$percent) {
        prcntii <- probObj$t
      } else if ((probObj$percent < percent) && (percent < probObj$percentP)) {
        h1 <- probObj$percentP - percent
        h2 <- percent - probObj$percent
        prcntii <- (h1 * probObj$t + h2 * probObj$tP) / (h1 + h2)
      } else if (percent == probObj$percentP) {
        prcntii <- probObj$tP
      }
      rnprcnt0[ii] <- as.numeric(prcntii)
    } else {
      probObj$isNull <- TRUE
    }
  }
  # smoothing output
  rnprcnt <- tsEvaNanRunningMean(rnprcnt0, windowSize)
  stdError <- sd(rnprcnt0 - rnprcnt)
  return(list(rnprcnt = rnprcnt, stdError = stdError))
}


#' Calculate the running mean trend of a time series
#'
#' This function calculates the running mean trend of a given time series using a specified time window.
#'
#' @param timeStamps A vector of time stamps corresponding to the observations in the series.
#' @param series A vector of numeric values representing the time series.
#' @param timeWindow The length of the time window (in the same time units as the time stamps) used for calculating the running mean trend.
#'
#' @return A list containing the running mean trend series and the number of observations used for each running mean calculation.
#'
#' @examples
#' timeStamps <- c(1, 2, 3, 4, 5)
#' series <- c(10, 15, 20, 25, 30)
#' timeWindow <- 2
#' result <- tsEvaRunningMeanTrend(timeStamps, series, timeWindow)
#' result$trendSeries
#' result$nRunMn
#'
#' @import stats
#' @export
tsEvaRunningMeanTrend <- function(timeStamps, series, timeWindow) {
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  nRunMn <- ceiling(timeWindow / dt)
  trendSeries1 <- tsEvaNanRunningMean(series, nRunMn)

  # further smoothing
  trendSeries <- tsEvaNanRunningMean(trendSeries1, ceiling(nRunMn / 2))
  return(list(trendSeries = trendSeries, nRunMn = nRunMn))
}

#' Detrend a Time Series
#'
#' This function detrends a time series by subtracting the trend component from the original series.
#'
#' @param timeStamps A vector of time stamps for the time series.
#' @param series The original time series.
#' @param timeWindow The size of the moving window used to calculate the trend.
#' @param percent The percentile value used to calculate the trend for extreme values. Default is NA.
#' @param fast A logical value indicating whether to print additional information. Default is FALSE.
#'
#' @import methods
#' @return An object of class "tsTrend" with the following components:
#' \describe{
#'   \item{originSeries}{The original time series.}
#'   \item{detrendSeries}{The detrended time series.}
#'   \item{trendSeries}{The trend component of the time series.}
#'   \item{nRunMn}{The number of data points in the moving window used to calculate the trend.}
#' }
#'
#' @examples
#' timeStamps <- c(1, 2, 3, 4, 5)
#' series <- c(10, 12, 15, 13, 11)
#' timeWindow <- 3
#' detrended <- tsEvaDetrendTimeSeries(timeStamps, series, timeWindow)
#' detrended$detrendSeries
#' @export
tsEvaDetrendTimeSeries <- function(timeStamps, series, timeWindow, percent = NA, fast = F) {
  extremeLowThreshold <- -Inf
  methods::setClass(
    Class = "tsTrend",
    representation(
      originSeries = "numeric",
      detrendSeries = "numeric",
      trendSeries = "numeric",
      nRunMn = "numeric"
    )
  )

  trendSeries <- tsEvaRunningMeanTrend(timeStamps, series, timeWindow)
  statSeries <- series
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  dtx <- dt
  if (tdim == "hours") dtx <- dt / 24
  if (fast == F) print(paste0("timestep: ", dt, tdim, " | time window: ", round(trendSeries$nRunMn * dtx / 365.25), " years"))
  statSeries[which(statSeries < extremeLowThreshold)] <- NA
  if (!is.na(percent)) {
    print(paste0("trend on the ", percent, " percentile"))
    trendSeries_perc <- tsEvaNanRunningPercentiles(timeStamps, series, trendSeries$nRunMn, percent)
    meanTrperc <- mean(trendSeries_perc$rnprcnt)
    # trend of the extremes
    trendSeriesxp <- trendSeries_perc$rnprcnt / meanTrperc * mean(trendSeries$trendSeries)
    trendSused <- trendSeriesxp
  } else {
    trendSused <- trendSeries$trendSeries
  }
  detrendSeries <- statSeries - trendSused
  return(new("tsTrend",
    originSeries = statSeries,
    detrendSeries = detrendSeries,
    trendSeries = trendSused,
    nRunMn = trendSeries$nRunMn
  ))
}


#' Estimate Average Seasonality
#'
#' This function estimates the average seasonality of a time series based on the given parameters.
#'
#' @param timeStamps The time stamps of the time series.
#' @param seasonalitySeries The series representing the seasonality.
#' @param timeWindow The time window used for averaging the seasonality.
#' @param peaks A logical value indicating whether to consider peaks in the seasonality.
#'
#' @return A list containing the estimated regime and the seasonality series.
#'   - \code{regime}: The estimated regime of the time series.
#'   - \code{Seasonality}: A data frame containing the average and varying seasonality series.
#'     - \code{averageSeasonalitySeries}: The average seasonality series.
#'     - \code{varyingSeasonalitySeries}: The varying seasonality series.
#'
#' @examples
#' timeStamps <- seq(as.Date("2022-01-01"), as.Date("2022-12-31"), by = "day")
#' seasonalitySeries <- sin(2 * pi * seq(0, 1, length.out = length(timeStamps)))
#' result <- tsEstimateAverageSeasonality(timeStamps, seasonalitySeries, timeWindow = 30, peaks = TRUE)
#' plot(result$regime, type = "l", xlab = "Day", ylab = "Regime", main = "Estimated Regime")
#' plot(result$Seasonality$averageSeasonalitySeries, type = "l", xlab = "Day", ylab = "Seasonality", main = "Average Seasonality")
#' plot(result$Seasonality$varyingSeasonalitySeries, type = "l", xlab = "Day", ylab = "Seasonality", main = "Varying Seasonality")
#'@importFrom pracma interp1
#' @export
tsEstimateAverageSeasonality <- function(timeStamps, seasonalitySeries, timeWindow, peaks = T) {
  avgYearLength <- 365.2425
  nMonthInYear <- 12
  timeWindow <- timeWindow
  # to get the number of days
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  avgYearLength <- avgYearLength / dt

  timstamp <- unique(as.Date(timeStamps))
  nYear <- round(length(timstamp) / avgYearLength)
  tw <- round(timeWindow / avgYearLength)
  avgMonthLength <- avgYearLength / nMonthInYear

  firstTmStmp <- timeStamps[1]
  lastTmStmp <- timeStamps[length(timeStamps)]
  mond <- lubridate::month(timeStamps)
  caca <- diff(mond)
  mony <- tsibble::yearmonth(timeStamps)
  monthTmStampStart <- c(timeStamps[1], timeStamps[which(diff(mond) != 0) + 1])
  monthTmStampEnd <- c(timeStamps[which(diff(mond) != 0)], timeStamps[length(timeStamps)])

  seasonalitySeries <- data.frame(time = timeStamps, series = seasonalitySeries, month = mony, mond = mond)

  if (peaks == T) {
    grpdSsn_ <- aggregate(seasonalitySeries$series,
      by = list(month = seasonalitySeries$month),
      FUN = function(x) sum(x, na.rm = T)
    )
  } else {
    grpdSsn_ <- aggregate(seasonalitySeries$series,
      by = list(month = seasonalitySeries$month),
      FUN = function(x) mean(x, na.rm = T)
    )
  }
  nYears <- ceiling(length(grpdSsn_$month) / nMonthInYear)
  grpdSsn <- matrix(NaN, nYears * nMonthInYear, 1)
  grpdSsn[1:length(grpdSsn_$x)] <- grpdSsn_$x
  grpdSsnMtx <- matrix(grpdSsn, nrow = 12, byrow = F)
  mnSsn_ <- rowMeans(grpdSsnMtx)

  # variations of the seasonality
  mnSsn_t <- c()
  for (y in (1 + tw / 2):(nYears - tw / 2)) {
    grpdSsnMti <- grpdSsnMtx[, c((y - tw / 2):(y + tw / 2))]
    twz <- round(timeWindow / 365 / dt)
    twv <- max(2, twz)
    mnSsn_i <- rowMeans(grpdSsnMti[, c(1:twv)])
    it <- (1:nMonthInYear)
    x <- it / 6. * pi
    dx <- pi / 6.
    a0 <- mean(mnSsn_)
    a1 <- 1 / pi * sum(cos(x) * mnSsn_i) * dx
    b1 <- 1 / pi * sum(sin(x) * mnSsn_i) * dx
    a2 <- 1 / pi * sum(cos(2 * x) * mnSsn_i) * dx
    b2 <- 1 / pi * sum(sin(2 * x) * mnSsn_i) * dx
    a3 <- 1 / pi * sum(cos(3 * x) * mnSsn_i) * dx
    b3 <- 1 / pi * sum(sin(3 * x) * mnSsn_i) * dx
    # mnSsni = a0  +  ( a1*cos(x) + b1*sin(x) )  +  ( a2*cos(2*x) + b2*sin(2*x) );
    mnSsni <- a0 + (a1 * cos(x) + b1 * sin(x)) + (a2 * cos(2 * x) + b2 * sin(2 * x)) + (a3 * cos(3 * x) + b3 * sin(3 * x))
    mnSsn_t <- c(mnSsn_t, mnSsni)
  }
  # estimating the first 2 fourier components
  it <- (1:nMonthInYear)
  x <- it / 6. * pi
  dx <- pi / 6.
  a0 <- mean(mnSsn_)
  a1 <- 1 / pi * sum(cos(x) * mnSsn_) * dx
  b1 <- 1 / pi * sum(sin(x) * mnSsn_) * dx
  a2 <- 1 / pi * sum(cos(2 * x) * mnSsn_) * dx
  b2 <- 1 / pi * sum(sin(2 * x) * mnSsn_) * dx
  a3 <- 1 / pi * sum(cos(3 * x) * mnSsn_) * dx
  b3 <- 1 / pi * sum(sin(3 * x) * mnSsn_) * dx
  mnSsn <- a0 + (a1 * cos(x) + b1 * sin(x)) + (a2 * cos(2 * x) + b2 * sin(2 * x)) + (a3 * cos(3 * x) + b3 * sin(3 * x))
  pt <- matrix(1, nYears, 1)
  monthAvgVec <- rep(mnSsn, nYears)
  imnth <- c(0:(length(monthAvgVec) - 1))
  avgTmStamp <- firstTmStmp + avgMonthLength / 2. + imnth * avgMonthLength
  imReg <- c(0:(length(mnSsn) - 1))
  regimeTmStamp <- 1 + avgMonthLength / 2 + imReg * avgMonthLength

  # adding first and last times
  monthAvgVec <- c(monthAvgVec[1], monthAvgVec, monthAvgVec[length(monthAvgVec)])
  avgTmStamp <- c(firstTmStmp, avgTmStamp, max(monthTmStampEnd))

  regimeVec <- c(mnSsn[1], mnSsn, mnSsn[length(mnSsn)])
  regimeTmStamp <- c(1, regimeTmStamp, 365)

  # for the varying seasonality
  # add the missiong 3 years at the end and beginning
  monthAvgVex <- c(rep(mnSsn_t[1:12], twz / 2), mnSsn_t, rep(mnSsn_t[(length(mnSsn_t) - 11):length(mnSsn_t)], twz / 2))
  monthAvgVex <- c(monthAvgVex[1], monthAvgVex, monthAvgVex[length(monthAvgVex)])
  avgTmStamp <- as.numeric(avgTmStamp)
  timeStampsN <- as.numeric(timeStamps)
  regime <- pracma::interp1(regimeTmStamp, regimeVec, c(1:365), method = "cubic")
  averageSeasonalitySeries <- pracma::interp1(avgTmStamp, monthAvgVec, timeStampsN, method = "cubic")
  varyingSeasonalitySeries <- pracma::interp1(avgTmStamp, monthAvgVex, timeStampsN, method = "cubic")
  return(list(regime = regime, Seasonality = data.frame(averageSeasonalitySeries = averageSeasonalitySeries, varyingSeasonalitySeries = varyingSeasonalitySeries)))
}


#' Calculate the return period of low flow based on a threshold and window size
#'
#' This function calculates the return period of low flow for a given time series
#' based on a threshold and window size. It uses a sliding window approach to
#' count the number of values below the threshold within each window, and then
#' calculates the return period based on the proportion of values below the
#' threshold. Assumes that the input data has a 7 days timestep.
#'
#' @param series The time series data.
#' @param threshold The threshold value for low flow.
#' @param windowSize The size of the sliding window.
#'
#' @return A data frame with two columns: "time" representing the time points
#'         corresponding to the sliding windows, and "RP" representing the
#'         calculated return period of low flow.
#'
#' @examples
#' series <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' threshold <- 5
#' windowSize <- 3
#' tsEvaNanRunnigBlowTh(series, threshold, windowSize)
#'
#' @export
tsEvaNanRunnigBlowTh <- function(series, threshold, windowSize) {
  halfWindow <- floor(windowSize / 2) / 28
  numDays <- seq(1, (length(series)))
  shortDay <- seq(1, length(series), by = 28)
  myd <- numeric(length(shortDay))
  ls <- numeric(length(shortDay))
  shortSerie <- series[shortDay]

  for (ii in seq(1, length(myd))) {
    startIdx <- max(1, (ii - 1) + 1 - halfWindow)
    endIdx <- min(length(shortSerie), ii + halfWindow)
    windowSubset <- shortSerie[startIdx:endIdx]
    myd[ii] <- sum(!is.na(windowSubset) & windowSubset <= threshold)
    ls[ii] <- length(windowSubset)
  }
  myd <- myd
  # smoothing with running mean
  mydsmooth <- tsEvaNanRunningMean(myd, halfWindow)
  p <- mydsmooth / ls
  # return period of low flow
  t <- 1 / (p * 52)
  return(data.frame(time = shortDay, RP = t))
}
