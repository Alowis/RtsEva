# Return Levels plots-------------

#' tsEvaPlotReturnLevelsGEVFromAnalysisObj
#'
#' \code{tsEvaPlotReturnLevelsGEVFromAnalysisObj} is a function that plots the
#' return levels for a Generalized Extreme Value (GEV) distribution using the
#' parameters obtained from an analysis object. It considers non-stationarity
#'  by considering time-varying parameters and their associated standard errors.
#'
#' @param nonStationaryEvaParams The non-stationary parameters obtained from the
#'  analysis object.
#' @param stationaryTransformData The stationary transformed data obtained from
#' the analysis object.
#' @param timeIndex The index at which the time-varying analysis should be
#' estimated.
#' @param trans The transformation used to fit the EVD. Can be "ori" for no
#' transformation or "rev" for reverse transformation.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return Plot 1: RLtstep: return level curve with confidence interval for the
#'  selected timeIndex
#'         Plot 2: beam: beam of return level curve for all with highlited curve
#'          for selected timeIndex
#'
#' @references
#' Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L., Besio,
#'  G., and Alfieri, L. (2016). The transformed-stationary approach: a generic
#'  and simplified methodology for non-stationary extreme value analysis.
#'   \emph{Hydrology and Earth System Sciences}, 20, 3527-3547.
#'    doi:10.5194/hess-20-3527-2016.
#'
#' @import ggplot2
#' @importFrom lubridate yday month
#' @importFrom texmex pgev
#' @docType methods
#' @name tsEvaPlotReturnLevelsGEVFromAnalysisObj
#' @export
tsEvaPlotReturnLevelsGEVFromAnalysisObj <- function(nonStationaryEvaParams,
                                                    stationaryTransformData,
                                                    timeIndex, trans, ...) {
  varargin <- NULL
  varargin <- list(...)
  args <- list(
    minReturnPeriodYears = 2,
    maxReturnPeriodYears = 1000,
    xlabel = "return period (years)",
    ylabel = "return levels (mm)",
    ylim = NULL,
    dtSampleYears = NULL, # one year
    ax = NULL
  )
  # Update args with passed in arguments
  varargin <- tsEasyParseNamedArgs(varargin, args)

  epsilon <- nonStationaryEvaParams[[1]]$parameters$epsilon
  sigma <- mean(nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex])
  mu <- mean(nonStationaryEvaParams[[1]]$parameters$mu[timeIndex])
  dtSampleYears <- nonStationaryEvaParams[[1]]$parameters$timeDeltaYears
  epsilonStdErr <- nonStationaryEvaParams[[1]]$paramErr$epsilonErr
  sigmaStdErr <- mean(nonStationaryEvaParams[[1]]$paramErr$sigmaErr[timeIndex])
  muStdErr <- mean(nonStationaryEvaParams[[1]]$paramErr$muErr[timeIndex])

  amax <- nonStationaryEvaParams[[1]]$parameters$annualMax
  monmax <- nonStationaryEvaParams[[1]]$parameters$monthlyMax
  amaxID <- nonStationaryEvaParams[[1]]$parameters$annualMaxIndx
  monmaxID <- nonStationaryEvaParams[[1]]$parameters$monthlyMaxIndx

  timeStamps <- stationaryTransformData$timeStamps
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  if (dt >= 1) {
    timeStamps <- stationaryTransformData$timeStamps
    trendAMax <- stationaryTransformData$trendSeries[amaxID]
    stdAMax <- stationaryTransformData$stdDevSeries[amaxID]
    amaxCor <- (amax - trendAMax) / stdAMax

    trendMax <- stationaryTransformData$trendSeries[monmaxID]
    stdMax <- stationaryTransformData$stdDevSeries[monmaxID]
    momaxCor <- (monmax - trendMax) / stdMax
  } else {
    timeStamps <- stationaryTransformData$timeStampsDay
    trendAMax <- stationaryTransformData$trendSeriesOr[amaxID]
    stdAMax <- stationaryTransformData$stdDevSeriesOr[amaxID]
    amaxCor <- (amax - trendAMax) / stdAMax

    trendMax <- stationaryTransformData$trendSeriesOr[monmaxID]
    stdMax <- stationaryTransformData$stdDevSeriesOr[monmaxID]
    momaxCor <- (monmax - trendMax) / stdMax
  }

  tstamps <- timeStamps[timeIndex]

  monmaxday <- lubridate::yday(stationaryTransformData$timeStamps[monmaxID])
  monmaxm <- lubridate::month(stationaryTransformData$timeStamps[monmaxID])
  xday <- lubridate::yday(stationaryTransformData$timeStamps[timeIndex]) * (2 * pi / 365.25)
  seasonalityTimeWindow <- 30.4
  maxD <- sqrt(2 - 2 * cos((2 * pi / 365.25) - seasonalityTimeWindow * (pi / 365.25)))
  seasonRatio <- seasonalityTimeWindow / 365.25
  theta <- monmaxday * (2 * pi / 365.25)
  D <- sqrt(2 - 2 * cos(xday - theta))
  sel <- which(month(stationaryTransformData$timeStamps[timeIndex]) == monmaxm)
  min(monmaxday[sel]) - lubridate::yday(stationaryTransformData$timeStamps[timeIndex])

  nYears <- length(amaxCor)
  rlvlamax <- empdis(amaxCor, nYears)
  rlvlamax$QNS <- amax[order(amax)]
  rlvlamax$Idt <- stationaryTransformData$timeStamps[amaxID][order(amax)]

  rlvlmmax <- empdis(momaxCor, nYears)
  rlvlmmax$QNS <- monmax[order(monmax)]

  rlvlmmax$gev.p <- (texmex::pgev(rlvlmmax$QNS, mu, sigma, epsilon, lower.tail = F))
  rlvlmmax$RP.gev <- 1 / (12 * rlvlmmax$gev.p)
  rlvmaxf <- rlvlmmax
  rlvmaxf$Idt <- stationaryTransformData$timeStamps[monmaxID][order(monmax)]

  if (nonStationaryEvaParams[[1]]$parameters$timeDeltaYears < 1) {
    rlvmax <- rlvmaxf
  } else {
    rlvmax <- rlvlamax
  }

  RLtstep <- tsEvaPlotReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, rlvmax, tstamps, trans, dtSampleYears)
  print(RLtstep)

  beam <- tsEvaPlotAllRLevelsGEV(nonStationaryEvaParams, stationaryTransformData, rlvmax, timeIndex, timeStamps, tstamps, trans, varargin)
  print(beam)
}

#' tsEvaPlotReturnLevelsGPDFromAnalysisObj
#'
#' \code{tsEvaPlotReturnLevelsGPDFromAnalysisObj} is a function that plots the return levels for a Generalized Pareto Distribution (GPD) using the parameters obtained from an analysis object. It considers non-stationarity by considering time-varying parameters and their associated standard errors.
#'
#' @param nonStationaryEvaParams The non-stationary parameters obtained from the analysis object.
#' @param stationaryTransformData The stationary transformed data obtained from the analysis object.
#' @param timeIndex The index at which the time-varying analysis should be estimated.
#' @param trans The transformation used to fit the EVD. Can be "ori" for no transformation or "rev" for reverse transformation.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return Plot 1: RLtstep: return level curve with confidence interval for the selected timeIndex
#'         Plot 2: beam: beam of return level curve for all with highlited curve for selected timeIndex
#'
#' @references
#' Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L., Besio, G., and Alfieri, L. (2016). The transformed-stationary approach: a generic and simplified methodology for non-stationary extreme value analysis. \emph{Hydrology and Earth System Sciences}, 20, 3527-3547. doi:10.5194/hess-20-3527-2016.
#'
#' @seealso [tsEvaPlotReturnLevelsGPD()] and [tsEvaPlotAllRLevelsGPD()]
#' @import ggplot2
#' @importFrom lubridate yday month
#' @name tsEvaPlotReturnLevelsGPDFromAnalysisObj
#' @aliases tsEvaPlotReturnLevelsGPDFromAnalysisObj
#' @export
tsEvaPlotReturnLevelsGPDFromAnalysisObj <- function(nonStationaryEvaParams,
                                                    stationaryTransformData,
                                                    timeIndex, trans, ope = F, ...) {
  varargin <- NULL
  varargin <- list(...)
  args <- list(
    minReturnPeriodYears = 0.5,
    maxReturnPeriodYears = 1000,
    xlabel = "return period (years)",
    ylabel = "return levels",
    ylim = NULL,
    dtSampleYears = NULL, # one year
    ax = NULL,
    trans = "rev"
  )
  # Update args with passed in arguments
  varargin <- tsEasyParseNamedArgs(varargin, args)

  epsilon <- nonStationaryEvaParams[[2]]$parameters$epsilon
  sigma <- mean(nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex])
  threshold <- mean(nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex])
  thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
  thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
  timeHorizonInYears <- round(as.numeric((thEnd - thStart) / 365.2425))
  nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks

  peax <- nonStationaryEvaParams[[2]]$parameters$peaks
  peaxID <- nonStationaryEvaParams[[2]]$parameters$peakID

  timeStamps <- stationaryTransformData$timeStamps
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  if (dt >= 1) {
    timeStamps <- stationaryTransformData$timeStamps
    trendPeaks <- stationaryTransformData$trendSeries[peaxID]
    stdPeaks <- stationaryTransformData$stdDevSeries[peaxID]
  } else {
    timeStamps <- stationaryTransformData$timeStampsDay
    trendPeaks <- stationaryTransformData$trendSeriesOr[peaxID]
    stdPeaks <- stationaryTransformData$stdDevSeriesOr[peaxID]
    dt <- 1
  }

  peaksCor <- (peax - trendPeaks) / stdPeaks
  tstamps <- timeStamps[timeIndex]


  epsilonStdErr <- nonStationaryEvaParams[[2]]$paramErr$epsilonErr
  sigmaStdErr <- mean(nonStationaryEvaParams[[2]]$paramErr$sigmaErr[timeIndex])
  thresholdStdErr <- mean(nonStationaryEvaParams[[2]]$paramErr$thresholdErr[timeIndex])

  peaxday <- lubridate::yday(stationaryTransformData$timeStamps[peaxID])
  xday <- lubridate::yday(stationaryTransformData$timeStamps[timeIndex]) * (2 * pi / 365.25)
  seasonalityTimeWindow <- 30.4 * 12
  maxD <- sqrt(2 - 2 * cos((2 * pi / 365.25) - seasonalityTimeWindow * (pi / 365.25)))
  seasonRatio <- seasonalityTimeWindow / 365.25
  theta <- peaxday * (2 * pi / 365.25)
  D <- sqrt(2 - 2 * cos(xday - theta))
  sel <- which(D < (maxD))

  nYears <- round(length(timeStamps) / 365.25 * dt)
  rlvlpeax <- empdis(peaksCor, nYears)
  rlvlpeax$QNS <- peax[order(peax)]
  rlvlpeax$Idt <- stationaryTransformData$timeStamps[peaxID][order(peax)]

  rlvlpeam <- empdis(peaksCor[sel], nYears * seasonRatio)
  rlvlpeam$QNS <- peax[sel][order(peax[sel])]
  rlvlpeam$Idt <- stationaryTransformData$timeStamps[peaxID[sel]][order(peax[sel])]


  if (nonStationaryEvaParams[[1]]$parameters$timeDeltaYears < 1) {
    rlvmax <- rlvlpeam
  } else {
    rlvmax <- rlvlpeax
  }

  RLtstep <- tsEvaPlotReturnLevelsGPD(epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr,
    nPeaks = nPeaks, timeHorizonInYears = timeHorizonInYears, rlvmax, tstamps, trans, varargin
  )
  print(RLtstep)

  beam <- tsEvaPlotAllRLevelsGPD(nonStationaryEvaParams, stationaryTransformData, rlvmax, timeIndex, timeStamps, tstamps, trans, varargin)
}

#' tsEvaPlotAllRLevelsGEV
#'
#' This function generates a beam plot of return levels for a Generalized Extreme Value (GEV) distribution based on the provided parameters and data.
#' The plot showcases the evolving relationship between return periods and return levels through time, allowing for visual analysis of extreme events and their probabilities.
#'
#' @param nonStationaryEvaParams A list of non-stationary evaluation parameters containing the GEV distribution parameters (epsilon, sigma, mu) and the time delta in years (dtSampleYears).
#' @param stationaryTransformData The stationary transformed data used for the analysis.
#' @param rlvmax The maximum return level data, including the return periods (haz.RP) and the actual return levels (QNS).
#' @param timeIndex The index of the time step used for analysis.
#' @param timeStamps The timestamps corresponding to the time steps in the analysis.
#' @param tstamps The timestamps used for labeling the plot.
#' @param trans The transformation used to fit the EVD, either "ori" (original) or "rev" (reverse).
#' @param ... Additional optional arguments for customizing the plot.
#'
#' @import ggplot2
#' @return A plot object showing the relationship between return periods and return levels for the GEV distribution at different timesteps.
#'
#' @examples
#' tsEvaPlotAllRLevelsGEV(nonStationaryEvaParams, stationaryTransformData, rlvmax, timeIndex, timeStamps, tstamps, mode = "inv", ylim = c(0, 100), ax = NULL)
#'
#' @seealso \code{\link{tsEvaComputeReturnLevelsGEV}}
#' @name tsEvaPlotAllRLevelsGEV
tsEvaPlotAllRLevelsGEV <- function(nonStationaryEvaParams, stationaryTransformData,
                                   rlvmax, timeIndex, timeStamps, tstamps,
                                   trans, varargin) {
  # Define default values for arguments
  epsilon <- nonStationaryEvaParams[[1]]$parameters$epsilon
  sigmav <- nonStationaryEvaParams[[1]]$parameters$sigma
  muv <- nonStationaryEvaParams[[1]]$parameters$mu
  dtSampleYears <- nonStationaryEvaParams[[1]]$parameters$timeDeltaYears
  timeStamps <- stationaryTransformData$timeStamps
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)

  args <- list(
    minReturnPeriodYears = 2,
    maxReturnPeriodYears = 1000,
    confidenceAreaColor = "lightgreen",
    confidenceBarColor = "darkgreen",
    returnLevelColor = "black",
    xlabel = "return period (years)",
    ylabel = "return levels",
    ylim = NULL,
    ax = NULL
  )

  # Update args with passed in arguments
  args <- tsEasyParseNamedArgs(varargin, args)
  minReturnPeriodYears <- args$minReturnPeriodYears
  maxReturnPeriodYears <- args$maxReturnPeriodYears

  # Compute return periods and levels
  returnPeriodsInYears <- 10^(seq(log10(minReturnPeriodYears), log10(maxReturnPeriodYears), by = 0.01))
  returnPeriodsInDts <- returnPeriodsInYears / dtSampleYears
  if (dtSampleYears < 1) {
    lts <- seq(1, length(sigmav), by = 32)
  }
  if (dtSampleYears >= 1) {
    lts <- seq(1, length(sigmav), by = 365.25 / dt)
  }
  rLevAll <- data.frame(matrix(ncol = 3, nrow = length(lts) * length(returnPeriodsInYears)))
  i <- 0
  lrp <- length(returnPeriodsInYears)
  for (ic in lts) {
    sigmax <- sigmav[ic]
    mux <- muv[ic]
    returnLevels <- tsEvaComputeReturnLevelsGEV(epsilon, sigmax, mux, epsilon, sigmax, mux, returnPeriodsInDts)$returnLevels
    if (trans == "inv") {
      returnLevels <- 1 / returnLevels
      rlvmax$Qreal <- 1 / rlvmax$QNS
    } else if (trans == "rev") {
      returnLevels <- -returnLevels
      rlvmax$Qreal <- -rlvmax$QNS
    } else if (trans == "lninv") {
      returnLevels <- 1 / exp(returnLevels)
      rlvmax$Qreal <- 1 / exp(rlvmax$QNS)
    } else {
      rlvmax$Qreal <- rlvmax$QNS
    }
    tic <- rep(as.Date(timeStamps[ic]), length(returnPeriodsInYears))
    Rper <- returnPeriodsInYears
    rLev <- data.frame(tic, Rper, as.vector(returnLevels))
    rLevAll[c((i * lrp + 1):(i * lrp + lrp)), ] <- rLev
    i <- i + 1
  }

  names(rLevAll) <- c("timeStamps", "returnPeriodsInYears", "returnLevels")
  rLevAll$timeStamps <- as.Date(rLevAll$timeStamps, origin = "1970-01-01")

  # Compute upper and lower bounds of return levels
  # Define y-axis limits
  if (!is.null(args$ylim)) {
    maxRL <- max(args$ylim)
    minRL <- min(args$ylim)
  } else if (trans == "inv") {
    maxRL <- round(max(rLevAll[, 3]) / 2) * 2 + 0.1
    minRL <- round(min(rLevAll[, 3]) / 2) * 1 - 0.1
  } else if (trans == "lninv") {
    maxRL <- round(max(rLevAll[, 3]) / 2) * 2 + 0.1
    minRL <- round(min(rLevAll[, 3]) / 2) * 1 - 0.1
  } else if (trans == "rev") {
    maxRL <- round(max(rLevAll[, 3]) * 100) / 100 + 0.01
    minRL <- round(min(rLevAll[, 3]) * 100) / 100 - 0.01
  } else {
    maxRL <- round(max(rLevAll[, 3]) / 10) * 10 + 5
    minRL <- round(min(rLevAll[, 3]) / 10) * 10 - 5
  }
  # Create plot
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))
  gpd.palette <- grDevices::colorRampPalette(c(
    "#F6FF33", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a",
    "#e31a1c", "#bd0026", "#800026", "#2C110B"
  ), interpolate = "linear", bias = 1)
  names(rLevAll) <- c("timeStamps", "returnPeriodsInYears", "returnLevels")
  rLevAll$timeStamps <- as.Date(rLevAll$timeStamps, origin = "1970-01-01")

  if (dtSampleYears < 1) {
    rLevAll$ts <- month(as.Date(rLevAll$timeStamps, origin = "1970-01-01"))
    rlvmax$cg <- month(as.Date(rlvmax$Idt, origin = "1970-01-01"))
    gpd.palette <- grDevices::colorRampPalette(c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"),
      interpolate = "linear", bias = 1
    )
    title <- "seasonal GEV curve beam"
    legendx <- "Months"
    spc <- "red"
  }
  if (dtSampleYears >= 1) {
    rLevAll$ts <- lubridate::year(rLevAll$timeStamps)
    rlvmax$cg <- lubridate::year(as.Date(rlvmax$Idt, origin = "1970-01-01"))
    gpd.palette <- grDevices::colorRampPalette(c(
      "#F6FF33", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a",
      "#e31a1c", "#bd0026", "#800026", "#2C110B"
    ), interpolate = "linear", bias = 1)
    title <- paste0("GEV curve beam - ", tstamps)
    legendx <- "Time"
    spc <- "darkblue"
  }
  rLevAll$group <- yearmonth(rLevAll$timeStamps, origin = "1970-01-01")

  IndexCurve <- tsEvaComputeReturnLevelsGEV(epsilon, sigmav[timeIndex], muv[timeIndex], epsilon, sigmav[timeIndex], muv[timeIndex], returnPeriodsInDts)$returnLevels
  IndexCurve <- as.vector(IndexCurve)

  # transformations
  if (trans == "inv") {
    IndexCurve <- 1 / IndexCurve
  } else if (trans == "rev") {
    IndexCurve <- tsEvaComputeReturnLevelsGEV(epsilon, -sigmav[timeIndex], -muv[timeIndex], epsilon, -sigmav[timeIndex], -muv[timeIndex], returnPeriodsInDts)$returnLevels
    IndexCurve <- as.vector(IndexCurve)
  } else if (trans == "lninv") {
    IndexCurve <- 1 / exp(IndexCurve)
  }
  IndexCurve <- data.frame(returnPeriodsInYears = returnPeriodsInYears, returnLevels = IndexCurve)

  f <- ggplot2::ggplot(rLevAll, aes(x = returnPeriodsInYears, y = returnLevels, color = ts, group = group)) +
    ggtitle(title) +
    geom_line(lwd = 1.5, alpha = 0.2) +
    geom_line(data = IndexCurve, aes(x = returnPeriodsInYears, y = returnLevels, group = 1), colour = spc, lwd = 1.5, alpha = 1) +
    geom_point(data = rlvmax, aes(x = haz.RP, y = Qreal, fill = cg, group = cg), pch = 21, size = 3, stroke = 1.5, color = "darkblue") +
    scale_colour_gradientn(guide = "colourbar", colours = gpd.palette(100), legendx) +
    scale_fill_gradientn(guide = "none", colours = gpd.palette(100), legendx) +
    annotate("label", x = , minReturnPeriodYears * 2, y = 0.9 * maxRL, label = paste0("\U03B5 = ", as.character(round(epsilon, 3))), size = 8) +
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks, args$xlabel, limits = c(minReturnPeriodYears, maxReturnPeriodYears)) +
    scale_y_continuous(
      n.breaks = 10, limit = c(minRL, maxRL), args$ylabel
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 24),
      panel.grid.minor.x = element_line(linetype = 2)
    )
  f
  return(f)
}

#' tsEvaPlotAllRLevelsGPD
#'
#' This function generates a plot of return levels for a Generalized Pareto
#' Distribution (GPD) based on the provided parameters and data.
#' The plot showcases the evolving relationship between return periods and
#' return levels, allowing for visual analysis of extreme events
#'  and their probabilities.
#'
#' @param nonStationaryEvaParams A list of non-stationary evaluation parameters
#'  containing the GPD distribution parameters (epsilon, sigma, threshold),
#'  time horizon start and end (thStart, thEnd), time horizon in years
#'  (timeHorizonInYears), and number of peaks (nPeaks).
#' @param stationaryTransformData The stationary transformed data used for
#'  the analysis.
#' @param rlvmax The maximum return level data, including the return periods
#'  (haz.RP) and the actual return levels (QNS).
#' @param timeIndex The index of the time step used for analysis.
#' @param timeStamps The timestamps corresponding to the time steps in
#' the analysis.
#' @param tstamps The timestamps used for labeling the plot.
#' @param trans The transformation used to fit the EVD, either
#' "ori" (original) or "rev" (reverse).
#' @param ... Additional optional arguments for customizing the plot.
#'
#' @import ggplot2
#' @return A plot object showing the relationship between return periods and
#' return levels for the GPD distribution at different timesteps.
#'
#'
#' @examples
#' tsEvaPlotAllRLevelsGPD(nonStationaryEvaParams, stationaryTransformData,
#' rlvmax, timeIndex, timeStamps, tstamps, mode = "inv", ylim = c(0, 100),
#'  ax = NULL)
#'

#' @seealso \code{\link{tsEvaComputeReturnLevelsGPD}}
#' @export
tsEvaPlotAllRLevelsGPD <- function(nonStationaryEvaParams, stationaryTransformData,
                                   rlvmax, timeIndex, timeStamps, tstamps,
                                   trans, varargin) {
  # Define default values for arguments
  epsilon <- nonStationaryEvaParams[[2]]$parameters$epsilon
  sigmao <- nonStationaryEvaParams[[2]]$parameters$sigma
  thresholdv <- nonStationaryEvaParams[[2]]$parameters$threshold
  thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
  thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
  timeHorizonInYears <- round(as.numeric((thEnd - thStart) / 365.2425))
  nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
  npy <- nPeaks / timeHorizonInYears
  timeStamps <- stationaryTransformData$timeStamps
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)


  args <- list(
    minReturnPeriodYears = 0.5,
    maxReturnPeriodYears = 1000,
    confidenceAreaColor = "lightgreen",
    confidenceBarColor = "darkgreen",
    returnLevelColor = "black",
    xlabel = "return period (years)",
    ylabel = "return levels (m3/s)",
    ylim = NULL,
    # dtSampleYears = dtSampleYears, # one year
    ax = NULL
  )

  # Update args with passed in arguments
  varargin <- list(...)
  args <- tsEasyParseNamedArgs(varargin, args)
  # print(args)
  minReturnPeriodYears <- args$minReturnPeriodYears
  maxReturnPeriodYears <- args$maxReturnPeriodYears

  # Compute return periods and levels
  returnPeriodsInYears <- 10^(seq(log10(minReturnPeriodYears), log10(maxReturnPeriodYears), by = 0.01))
  if (npy > 5) {
    lts <- seq(1, length(sigmao), by = 32)
  }
  if (npy <= 5) {
    lts <- seq(1, length(sigmao), by = 365.25 / dt)
  }
  rLevAll <- data.frame(matrix(ncol = 3, nrow = length(lts) * length(returnPeriodsInYears)))
  i <- 0
  lrp <- length(returnPeriodsInYears)
  for (ic in lts) {
    sigmax <- sigmao[ic]
    thresholdx <- thresholdv[ic]
    returnLevels <- tsEvaComputeReturnLevelsGPD(epsilon, sigmax, thresholdx, epsilon, sigmax, thresholdx, nPeaks = nPeaks, sampleTimeHorizon = timeHorizonInYears, returnPeriodsInYears)$returnLevels
    if (trans == "inv") {
      returnLevels <- 1 / returnLevels
      rlvmax$Qreal <- 1 / rlvmax$QNS
    } else if (trans == "rev") {
      returnLevels <- -returnLevels
      rlvmax$Qreal <- -rlvmax$QNS
    } else if (trans == "lninv") {
      returnLevels <- 1 / exp(returnLevels)
      rlvmax$Qreal <- 1 / exp(rlvmax$QNS)
    } else {
      rlvmax$Qreal <- rlvmax$QNS
    }
    tic <- rep(as.Date(timeStamps[ic]), length(returnPeriodsInYears))
    Rper <- returnPeriodsInYears
    rLev <- data.frame(tic, Rper, returnLevels)
    rLevAll[c((i * lrp + 1):(i * lrp + lrp)), ] <- rLev
    i <- i + 1
  }

  names(rLevAll) <- c("timeStamps", "returnPeriodsInYears", "returnLevels")
  rLevAll$returnLevels[which(rLevAll$returnLevels < 0)] <- 0

  rLevAll$timeStamps <- as.Date(rLevAll$timeStamps, origin = "1970-01-01")

  # Compute upper and lower bounds of return levels
  # Define y-axis limits
  if (!is.null(args$ylim)) {
    maxRL <- max(args$ylim)
    minRL <- min(args$ylim)
  } else if (trans == "inv") {
    maxRL <- round(max(rLevAll[, 3]) / 2) * 2 + 0.1
    minRL <- round(min(rLevAll[, 3]) / 2) * 1 - 0.1
  } else if (trans == "rev") {
    maxRL <- round(max(rLevAll[, 3]) * 100) / 100 + 0.01
    minRL <- round(min(rLevAll[, 3]) * 100) / 100 - 0.01
  } else if (trans == "lninv") {
    maxRL <- round(max(rLevAll[, 3]) / 2) * 2 + 0.1
    minRL <- round(min(rLevAll[, 3]) / 2) * 1 - 0.1
  } else {
    maxRL <- round(max(rLevAll[, 3]) / 10) * 10 + 5
    minRL <- round(min(rLevAll[, 3]) / 10) * 10 - 5
  }

  # Create plot
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))
  if (npy > 6) {
    rLevAll$ts <- month(as.Date(rLevAll$timeStamps, origin = "1970-01-01"))
    rlvmax$cg <- month(as.Date(rlvmax$Idt, origin = "1970-01-01"))
    gpd.palette <- colorRampPalette(c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"),
      interpolate = "linear", bias = 1
    )
    title <- "seasonal GPD curve beam"
    legendx <- "Months"
  }
  if (npy <= 6) {
    rLevAll$ti <- lubridate::year(rLevAll$timeStamps)
    rlvmax$cg <- lubridate::year(as.Date(rlvmax$Idt, origin = "1970-01-01"))
    gpd.palette <- colorRampPalette(c(
      "#F6FF33", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a",
      "#e31a1c", "#bd0026", "#800026", "#2C110B"
    ), interpolate = "linear", bias = 1)
    title <- paste0("GPD curve beam - ", tstamps)
    legendx <- "Time"
  }

  rLevAll$group <- as.numeric(yearmonth(rLevAll$timeStamps))
  IndexCurve <- tsEvaComputeReturnLevelsGPD(epsilon, sigmao[timeIndex], thresholdv[timeIndex], epsilon, sigmao[timeIndex], thresholdv[timeIndex], nPeaks = nPeaks, sampleTimeHorizon = timeHorizonInYears, returnPeriodsInYears)$returnLevels
  if (trans == "inv") {
    IndexCurve <- 1 / IndexCurve
  } else if (trans == "rev") {
    IndexCurve <- -IndexCurve
  } else if (trans == "lninv") {
    IndexCurve <- 1 / exp(IndexCurve)
  }
  IndexCurve <- data.frame(returnPeriodsInYears = returnPeriodsInYears, returnLevels = IndexCurve)
  IndexCurve$returnLevels[which(IndexCurve$returnLevels < 0)] <- 0

  f <- ggplot(rLevAll, aes(x = returnPeriodsInYears, y = returnLevels, color = ti, group = group)) +
    ggtitle(title) +
    geom_line(lwd = 1.5, alpha = 0.5) +
    geom_line(data = IndexCurve, aes(x = returnPeriodsInYears, y = returnLevels, group = 1), colour = "darkblue", lwd = 1.5, alpha = 1) +
    geom_point(data = rlvmax, aes(x = haz.RP, y = Qreal, fill = cg, group = cg), pch = 21, size = 3, stroke = 1.5, color = "darkblue") +
    scale_colour_gradientn(guide = "colourbar", colours = gpd.palette(100), legendx) +
    scale_fill_gradientn(guide = "none", colours = gpd.palette(100), "Time") +
    annotate("label", x = , maxReturnPeriodYears / 2, y = maxRL, label = paste0("\U03B5 = ", as.character(round(epsilon, 3))), size = 8) +
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks, args$xlabel, limits = c(minReturnPeriodYears, maxReturnPeriodYears)) +
    scale_y_continuous(
      n.breaks = 10, limit = c(minRL, maxRL), args$ylabel
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 24),
      panel.grid.minor.x = element_line(linetype = 2)
    )
  f
  return(f)
}


#' tsEvaPlotReturnLevelsGEV Function
#'
#' This function plots the return levels using the Generalized Extreme Value (GEV) distribution.
#'
#' @param epsilon The shape parameter of the GEV distribution.
#' @param sigma The scale parameter of the GEV distribution.
#' @param mu The location parameter of the GEV distribution.
#' @param epsilonStdErr The standard error of the shape parameter.
#' @param sigmaStdErr The standard error of the scale parameter.
#' @param muStdErr The standard error of the location parameter.
#' @param rlvmax A data frame containing the return levels of annual maxima.
#' @param tstamps The title for the plot.
#' @param trans The transformation used to fit the EVD, either "ori" (original) or "rev" (reverse).
#' @param ... Additional arguments to be passed to the function.
#'
#' @import ggplot2
#' @return A ggplot object representing the plot of return levels.
#' @seealso \code{\link{tsEvaComputeReturnLevelsGEV}} \code{\link{tsEvaPlorReturnLevelsGEVFromAnalysisObj}}
#' @examples
#' tsEvaPlotReturnLevelsGEV(
#'   epsilon = 0.1, sigma = 1, mu = 0, epsilonStdErr = 0.05, sigmaStdErr = 0.1, muStdErr = 0.05,
#'   rlvmax = data.frame(QNS = c(1, 2, 3), Qreal = c(10, 20, 30)),
#'   tstamps = "Return Levels", trans = "inv")
#'   @export
tsEvaPlotReturnLevelsGEV <- function(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr,
                                     muStdErr, rlvmax, tstamps, trans, ...) {
  varargin <- NULL
  varagin <- list(...)
  # Define default values for arguments
  args <- list(
    minReturnPeriodYears = 2,
    maxReturnPeriodYears = 1000,
    confidenceAreaColor = "lightgreen",
    confidenceBarColor = "darkgreen",
    returnLevelColor = "black",
    xlabel = "return period (years)",
    ylabel = "return levels",
    ylim = NULL,
    dtSampleYears = 1, # one year
    ax = NULL
  )

  # Update args with passed in arguments
  args <- tsEasyParseNamedArgs(varargin, args)
  minReturnPeriodYears <- args$minReturnPeriodYears
  maxReturnPeriodYears <- args$maxReturnPeriodYears

  # Compute return periods and levels
  returnPeriodsInYears <- 10^(seq(log10(minReturnPeriodYears), log10(maxReturnPeriodYears), by = 0.01))
  returnPeriodsInDts <- returnPeriodsInYears / args$dtSampleYears
  returnLevels <- tsEvaComputeReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriodsInDts)

  returnPeriodsInYears <- returnPeriodsInYears[which(!is.infinite(returnLevels$returnLevels))]
  returnPeriodsInDts <- returnPeriodsInDts[which(!is.infinite(returnLevels$returnLevels))]
  returnLevels$returnLevelsErr <- returnLevels$returnLevelsErr[which(!is.infinite(returnLevels$returnLevels))]
  returnLevels$returnLevels <- returnLevels$returnLevels[which(!is.infinite(returnLevels$returnLevels))]

  # Compute upper and lower bounds of return levels
  supRLCI <- returnLevels$returnLevels + returnLevels$returnLevelsErr
  infRLCI <- returnLevels$returnLevels - returnLevels$returnLevelsErr
  if (trans == "inv") {
    returnLevels$returnLevels <- 1 / returnLevels$returnLevels
    supRLCI <- 1 / supRLCI
    infRLCI <- 1 / infRLCI
    rlvmax$Qreal <- 1 / rlvmax$QNS
  } else if (trans == "rev") {
    returnLevels$returnLevels <- -returnLevels$returnLevels
    supRLCI <- -supRLCI
    infRLCI <- -infRLCI
    rlvmax$Qreal <- -rlvmax$QNS
  } else if (trans == "lninv") {
    returnLevels$returnLevels <- 1 / exp(returnLevels$returnLevels)
    supRLCI <- 1 / exp(supRLCI)
    infRLCI <- 1 / exp(infRLCI)
    rlvmax$Qreal <- 1 / exp(rlvmax$QNS)
  } else {
    rlvmax$Qreal <- rlvmax$QNS
  }

  # Define y-axis limits
  if (!is.null(args$ylim)) {
    maxRL <- max(args$ylim)
    minRL <- min(args$ylim)
  } else if (trans == "inv") {
    maxRL <- round(max(infRLCI) / 2) * 2 + 0.5
    minRL <- round(min(supRLCI) / 2) * 1 - 0.5
  } else if (trans == "rev") {
    maxRL <- round(max(infRLCI) * 100) / 100 + 0.01
    minRL <- round(min(supRLCI) * 100) / 100 - 0.01
  } else if (trans == "lninv") {
    maxRL <- round(max(infRLCI) / 2) * 2 + 0.5
    minRL <- round(min(supRLCI) / 2) * 1 - 0.5
  } else {
    maxRL <- round(max(supRLCI) / 10) * 10 + 5
    minRL <- round(min(supRLCI) / 10) * 10 - 5
  }

  # Create plot
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))

  gpd.palette <- colorRampPalette(c(
    "#F6FF33", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a",
    "#e31a1c", "#bd0026", "#800026", "#2C110B"
  ), interpolate = "linear", bias = 1)

  df <- data.frame(
    returnPeriodsInYears = returnPeriodsInYears, returnLevels = returnLevels$returnLevels,
    infRLCI = infRLCI, supRLCI = supRLCI
  )


  f <- ggplot(df, aes(x = returnPeriodsInYears, y = returnLevels)) +
    geom_ribbon(aes(ymin = infRLCI, ymax = supRLCI), fill = args$confidenceAreaColor, alpha = 0.5) +
    geom_line(color = args$returnLevelColor, lwd = 1.5) +
    geom_line(aes(y = supRLCI), color = args$confidenceBarColor, lwd = 1.2) +
    geom_line(aes(y = infRLCI), color = args$confidenceBarColor, lwd = 1.2) +
    geom_point(data = rlvmax, aes(x = haz.RP, y = Qreal, fill = lubridate::year(Idt)), pch = 21, size = 3, stroke = 1.5, color = "darkblue") +
    scale_fill_gradientn(guide = "colourbar", colours = gpd.palette(100), "Time") +
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks, args$xlabel, limits = c(minReturnPeriodYears, maxReturnPeriodYears)) +
    scale_y_continuous(
      n.breaks = 10, limit = c(minRL, maxRL), args$ylabel
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 24),
      panel.grid.minor.x = element_line(linetype = 2)
    ) +
    ggtitle(tstamps)
  f

  return(f)
}

#' Plot Return Levels using Generalized Pareto Distribution (GPD)
#'
#' This function plots the return levels using the Generalized Pareto Distribution (GPD).
#'
#' @param epsilon The shape parameter of the GPD.
#' @param sigma The scale parameter of the GPD.
#' @param threshold The threshold parameter of the GPD.
#' @param epsilonStdErr The standard error of the shape parameter.
#' @param sigmaStdErr The standard error of the scale parameter.
#' @param thresholdStdErr The standard error of the threshold parameter.
#' @param nPeaks The number of peaks used in the GPD estimation.
#' @param timeHorizonInYears The time horizon in years for the GPD estimation.
#' @param rlvmax A data frame containing the return levels of annual maxima.
#' @param tstamps The title for the plot.
#' @param trans The transformation type for the return levels.
#' @param ... Additional arguments to be passed to the function.
#'
#' @import ggplot2
#' @seealso \code{\link{tsEvaComputeReturnLevelsGPD}} \code{\link{tsEvaPlorReturnLevelsGPDFromAnalysisObj}}
#' @return A ggplot object representing the plot of return levels.
#'
#' @examples
#' tsEvaPlotReturnLevelsGPD(
#'   epsilon = 0.1, sigma = 1, threshold = 0, epsilonStdErr = 0.05, sigmaStdErr = 0.1,
#'   thresholdStdErr = 0.05, nPeaks = 10, timeHorizonInYears = 50,
#'   rlvmax = data.frame(QNS = c(1, 2, 3), Qreal = c(10, 20, 30)),
#'   tstamps = "Return Levels", trans = "inv")
#'   @export
tsEvaPlotReturnLevelsGPD <- function(epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, timeHorizonInYears, rlvmax, tstamps, trans, ...) {
  varargin <- list(...)

  # Define default values for arguments
  args <- list(
    minReturnPeriodYears = 0.5,
    maxReturnPeriodYears = 1000,
    confidenceAreaColor = "lightgreen",
    confidenceBarColor = "darkgreen",
    returnLevelColor = "black",
    xlabel = "return period (years)",
    ylabel = "return levels (mm)",
    ylim = NULL,
    dtSampleYears = 1, # one year
    ax = NULL
  )
  # Update args with passed in arguments
  # print(varargin)
  args <- tsEasyParseNamedArgs(varargin, args)
  # print(args)
  minReturnPeriodYears <- args$minReturnPeriodYears
  maxReturnPeriodYears <- args$maxReturnPeriodYears

  # Compute return periods and levels
  returnPeriods <- 10^(seq(log10(minReturnPeriodYears), log10(maxReturnPeriodYears), by = 0.01))


  returnLevels <- tsEvaComputeReturnLevelsGPD(epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr,
    nPeaks = nPeaks, sampleTimeHorizon = timeHorizonInYears, returnPeriods = returnPeriods
  )


  # retrunLevelsErrs <- tsEvaComputeReturnLevelsGEV(epsilonStdErr = epsilonStdErr, sigmaStdErr = sigmaStdErr, muStdErr = muStdErr, returnPeriodsInDts = returnPeriodsInDts)

  # Compute upper and lower bounds of return levels
  supRLCI <- returnLevels$returnLevels + returnLevels$returnLevelsErr
  infRLCI <- returnLevels$returnLevels - returnLevels$returnLevelsErr
  if (trans == "inv") {
    returnLevels$returnLevels <- 1 / returnLevels$returnLevels
    supRLCI <- 1 / supRLCI
    infRLCI <- 1 / infRLCI
    rlvmax$Qreal <- 1 / rlvmax$QNS
  } else if (trans == "rev") {
    returnLevels$returnLevels <- -returnLevels$returnLevels
    supRLCI <- -supRLCI
    infRLCI <- -infRLCI
    rlvmax$Qreal <- -rlvmax$QNS
  } else if (trans == "lninv") {
    returnLevels$returnLevels <- 1 / exp(returnLevels$returnLevels)
    supRLCI <- 1 / exp(supRLCI)
    infRLCI <- 1 / exp(infRLCI)
    rlvmax$Qreal <- 1 / exp(rlvmax$QNS)
  } else {
    rlvmax$Qreal <- rlvmax$QNS
  }

  # Define y-axis limits
  if (!is.null(args$ylim)) {
    maxRL <- round(max(args$ylim))
    minRL <- round(min(args$ylim))
  } else if (trans == "inv") {
    maxRL <- round(max(infRLCI) / 2) * 2 + 0.01
    minRL <- round(min(supRLCI) / 2) * 2 - 0.01
  } else if (trans == "rev") {
    maxRL <- round(max(infRLCI) * 100) / 100 + 0.01
    minRL <- round(min(supRLCI) * 100) / 100 - 0.01
  } else if (trans == "lminv") {
    maxRL <- round(max(infRLCI) / 2) * 2 + 0.01
    minRL <- round(min(supRLCI) / 2) * 2 - 0.01
  } else {
    maxRL <- round(max(supRLCI) / 5) * 5 + 5
    minRL <- round(min(infRLCI) / 5) * 5 - 5
  }
  gpd.palette <- colorRampPalette(c(
    "#F6FF33", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a",
    "#e31a1c", "#bd0026", "#800026", "#2C110B"
  ), interpolate = "linear", bias = 1)
  # remove 0 data from CI

  infRLCI[which(infRLCI < 0)] <- 0
  supRLCI[which(supRLCI < 0)] <- 0

  # Create plot
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))

  dfr <- data.frame(
    returnPeriods = returnPeriods, returnLevels = returnLevels$returnLevels,
    infRLCI = infRLCI, supRLCI = supRLCI
  )
  neg <- which(dfr$returnLevels < 0)
  if (length(neg) > 1) idm <- neg[1]
  dfr$returnLevels[which(dfr$returnLevels < 0)] <- 0
  f <- ggplot(dfr, mapping = aes(x = returnPeriods, y = returnLevels)) +
    geom_ribbon(aes(ymin = infRLCI, ymax = supRLCI), fill = args$confidenceAreaColor, alpha = 0.5) +
    geom_line(color = args$returnLevelColor, lwd = 1.5) +
    geom_line(aes(y = supRLCI), color = args$confidenceBarColor, lwd = 1.2) +
    geom_line(aes(y = infRLCI), color = args$confidenceBarColor, lwd = 1.2) +
    geom_point(data = rlvmax, aes(x = haz.RP, y = Qreal, fill = lubridate::year(Idt)), pch = 21, size = 3, stroke = 1.5, color = "darkblue") +
    scale_fill_gradientn(guide = "colourbar", colours = gpd.palette(100), "Years") +
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks, args$xlabel, limits = c(minReturnPeriodYears, maxReturnPeriodYears)) +
    scale_y_continuous(
      n.breaks = 10, limit = c(minRL, maxRL), args$ylabel
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 24),
      panel.grid.minor.x = element_line(linetype = 2)
    ) +
    ggtitle(tstamps)
  if (length(neg) > 1) {
    f <- f + geom_vline(xintercept = dfr$returnPeriods[idm], col = "red", lwd = 2)
  }

  f
  # find a better way to save the plot
  # if(!is.null(f))
  #   ggsave(f, plot = last_plot(), device = "pdf")

  return(f)
  # consider adding the observations to that plot
}



#' Generate GEV Image Scatter Plot from Analysis Object
#'
#' This function generates a GEV (Generalized Extreme Value) image scatter plot from an analysis object.
#'
#' @param Y The input data.
#' @param nonStationaryEvaParams A list of non-stationary evaluation parameters.
#' @param stationaryTransformData The stationary transform data.
#' @param trans The transformation method.
#' @param ... Additional arguments.
#'
#' @import ggplot2
#' @return The GEV image scatter plot.
#'
#' @examples
#' # Example usage of tsEvaPlotGEVImageScFromAnalysisObj function
#' tsEvaPlotGEVImageScFromAnalysisObj(Y, nonStationaryEvaParams, stationaryTransformData, trans)
#'
#' @export
tsEvaPlotGEVImageScFromAnalysisObj <- function(Y, nonStationaryEvaParams, stationaryTransformData, trans, ...) {
  varargin <- NULL
  varargin <- list(...)

  timeStamps <- stationaryTransformData$timeStamps
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  if (dt == 1) {
    timeStamps <- stationaryTransformData$timeStamps
  } else {
    timeStamps <- stationaryTransformData$timeStampsDay
  }

  epsilon <- nonStationaryEvaParams[[1]]$parameters$epsilon
  serix <- stationaryTransformData$nonStatSeries
  sigma <- nonStationaryEvaParams[[1]]$parameters$sigma
  mu <- nonStationaryEvaParams[[1]]$parameters$mu
  dtSampleYears <- nonStationaryEvaParams[[1]]$parameters$timeDeltaYears
  returnPeriodsInDts <- 1 / dtSampleYears

  amax <- nonStationaryEvaParams[[1]]$parameters$annualMax
  monmax <- nonStationaryEvaParams[[1]]$parameters$monthlyMax
  amaxID <- nonStationaryEvaParams[[1]]$parameters$annualMaxIndx
  monmaxID <- nonStationaryEvaParams[[1]]$parameters$monthlyMaxIndx

  tamax <- stationaryTransformData$timeStamps[amaxID]
  tmonmax <- stationaryTransformData$timeStamps[monmaxID]

  amaxplot <- data.frame(time = tamax, value = amax)
  monmaxplot <- data.frame(time = tmonmax, value = monmax)


  if (nonStationaryEvaParams[[1]]$parameters$timeDeltaYears <= 1) {
    maxObs <- amaxplot
  } else {
    maxObs <- monmaxplot
  }

  plotbg <- tsEvaPlotGEVImageSc(Y, timeStamps, serix, epsilon, sigma, mu, returnPeriodsInDts, maxObs, mode, varargin)
  print(plotbg)
  return(plotbg)
}


#' Plot GPD Image Score from Analysis Object
#'
#' This function plots the GPD (Generalized Pareto Distribution) image score from an analysis object.
#'
#' @param Y The input data.
#' @param nonStationaryEvaParams A list containing non-stationary evaluation parameters.
#' @param stationaryTransformData A data frame containing stationary transform data.
#' @param trans The transformation method to be applied to the data.
#' @param ... Additional arguments to be passed to the \code{\link{tsEvaPlotGPDImageSc}} function.
#'
#' @import ggplot2
#' @return The plot object.
#'
#' @details This function takes the input data \code{Y}, non-stationary evaluation parameters \code{nonStationaryEvaParams},
#' stationary transform data \code{stationaryTransformData}, transformation method \code{trans}, and additional arguments \code{...}.
#' It then updates the arguments with the passed-in values, calculates the time stamps, and performs necessary transformations.
#' Finally, it plots the GPD image score using the \code{\link{tsEvaPlotGPDImageSc}} function and returns the plot object.
#'
#' @seealso \code{\link{tsEvaPlotGPDImageSc}}
#' @seealso \code{\link{stationaryTransformData}}
#' @seealso \code{\link{nonStationaryEvaParams}}
#'
#' @examples
#' # Example usage of tsEvaPlotGPDImageScFromAnalysisObj
#' tsEvaPlotGPDImageScFromAnalysisObj(Y, nonStationaryEvaParams, stationaryTransformData, trans)
#'
#' @export
tsEvaPlotGPDImageScFromAnalysisObj <- function(Y, nonStationaryEvaParams, stationaryTransformData, trans, ...) {
  varargin <- NULL
  varargin <- list(...)

  # Update args with passed in arguments
  timeStamps <- stationaryTransformData$timeStamps
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  if (dt >= 1) {
    timeStamps <- stationaryTransformData$timeStamps
  } else {
    timeStamps <- stationaryTransformData$timeStampsDay
  }
  serix <- stationaryTransformData$nonStatSeries
  epsilon <- nonStationaryEvaParams[[2]]$parameters$epsilon
  sigma <- nonStationaryEvaParams[[2]]$parameters$sigma
  threshold <- nonStationaryEvaParams[[2]]$parameters$threshold
  Y <- Y

  peax <- nonStationaryEvaParams[[2]]$parameters$peaks
  peaxID <- nonStationaryEvaParams[[2]]$parameters$peakID

  peaktime <- stationaryTransformData$timeStamps[peaxID]
  peakplot <- data.frame(time = peaktime, value = peax)
  plotbg <- tsEvaPlotGPDImageSc(Y, as.Date(timeStamps), serix, epsilon, sigma, threshold, peakplot, trans, varargin)
  print(plotbg)
  return(plotbg)
}



#' tsEvaPlotGPDImageSc function
#'
#' This function generates a plot of the Generalized Pareto Distribution (GPD) using the provided data.
#'
#' @param Y A vector of values.
#' @param timeStamps A vector of timestamps corresponding to the values.
#' @param serix A vector of series values.
#' @param epsilon A numeric value representing the shape parameter of the GPD.
#' @param sigma A vector of standard deviation values.
#' @param threshold A vector of threshold values.
#' @param peakplot A data frame containing peak values and their corresponding timestamps.
#' @param trans A character string indicating the transformation to be applied to the data.
#' @param varargin Additional optional arguments.
#'
#' @import ggplot2 scales
#' @return A ggplot object representing the GPD plot.
#'
#' @examples
#' # Example usage of tsEvaPlotGPDImageSc function
#' plot <- tsEvaPlotGPDImageSc(Y, timeStamps, serix, epsilon, sigma, threshold, peakplot, trans, ...)
#' print(plot)
#'
#' @export
tsEvaPlotGPDImageSc <- function(Y, timeStamps, serix, epsilon, sigma, threshold, peakplot, trans, varargin) {
  avgYearLength <- 365.2425
  nyears <- as.numeric(round((max(timeStamps) - min(timeStamps)) / avgYearLength))
  nelmPerYear <- length(timeStamps) / nyears

  args <- list(
    nPlottedTimesByYear = min(12, round(nelmPerYear)),
    ylabel = "Values",
    xlabel = "Date",
    minYear = 1950,
    maxYear = 2021,
    axisFontSize = 20,
    labelFontSize = 22,
    Title = ""
  )
  args <- tsEasyParseNamedArgs(varargin, args)

  minTS <- as.Date(paste(args$minYear, "-01-01", sep = ""))
  maxTS <- as.Date(paste(args$maxYear, "-01-01", sep = ""))
  sigma <- sigma[(timeStamps >= minTS) & (timeStamps <= maxTS)]
  threshold <- threshold[(timeStamps >= minTS) & (timeStamps <= maxTS)]
  timeStamps <- timeStamps[(timeStamps >= minTS) & (timeStamps <= maxTS)]

  serie <- data.frame(timeStamps, serix)

  L <- length(timeStamps)
  minTS <- as.Date(timeStamps[1])
  maxTS <- as.Date(timeStamps[length(timeStamps)])
  npdf <- as.numeric(ceiling(((maxTS - minTS) / avgYearLength) * args$nPlottedTimesByYear))
  navg <- 1
  # navg <- floor(L/npdf)

  plotSLength <- npdf * navg
  timeStamps_plot <- seq(minTS, maxTS, length.out = plotSLength)

  if (length(epsilon) == 1) {
    epsilon0 <- rep(epsilon, npdf)
  } else {
    epsilon_ <- rep(NA, npdf * navg)
    epsilon_[1:L] <- epsilon
    epsilonMtx <- matrix(epsilon_, navg, ncol = npdf)
    epsilon0 <- apply(epsilonMtx, 2, mean, na.rm = TRUE)
  }

  sigma_ <- approx(timeStamps, sigma, timeStamps_plot)$y
  sigmaMtx <- matrix(sigma_, navg, ncol = npdf)
  if (length(sigmaMtx) > 1) {
    sigma0 <- apply(sigmaMtx, 2, mean, na.rm = TRUE)
  } else {
    sigma0 <- sigmaMtx
  }

  threshold_ <- approx(timeStamps, threshold, timeStamps_plot)$y
  thresholdMtx <- matrix(threshold_, navg, ncol = npdf)
  if (length(thresholdMtx) > 1) {
    threshold0 <- apply(thresholdMtx, 2, mean, na.rm = TRUE)
  } else {
    threshold0 <- thresholdMtx
  }


  # extremesRange = c(9,100)
  # Y <- c(seq(min(extremesRange),max(extremesRange),length.out=300))
  epsilon0 <- as.numeric(epsilon0)
  sigma0 <- as.numeric(sigma0)
  threshold0 <- as.numeric(threshold0)
  gridEps <- expand.grid(Y, epsilon0)
  gridSig <- expand.grid(Y, sigma0)
  gridthreshold <- expand.grid(Y, threshold0)
  gridTime <- expand.grid(Y, timeStamps_plot)


  # colnames(grid) <- c("Y", "epsilon", "sigma", "threshold")
  if (trans == "inv") {
    Ybis <- seq(min(1 / Y), max(1 / Y), length.out = length(Y))
    pdf <- texmex::dgpd(1 / Ybis, gridSig$Var2, gridEps$Var2, gridthreshold$Var2)
    datap <- data.frame(
      timeStamps = as.Date(gridTime$Var2), extremeValues = Ybis,
      pdf = pdf
    )
  } else if (trans == "rev") {
    Ybis <- seq(min(-Y), max(-Y), length.out = length(Y))
    pdf <- texmex::dgpd(-Ybis, gridSig$Var2, gridEps$Var2, gridthreshold$Var2)
    datap <- data.frame(
      timeStamps = as.Date(gridTime$Var2), extremeValues = Ybis,
      pdf = pdf
    )
  } else if (trans == "lninv") {
    Ybis <- seq(min(1 / exp(Y)), max(1 / exp(Y)), length.out = length(Y))
    pdf <- texmex::dgpd(-log(Ybis), gridSig$Var2, gridEps$Var2, gridthreshold$Var2)
    datap <- data.frame(
      timeStamps = as.Date(gridTime$Var2), extremeValues = Ybis,
      pdf = pdf
    )
  } else {
    pdf <- texmex::dgpd(Y, gridSig$Var2, gridEps$Var2, gridthreshold$Var2)
    datap <- data.frame(
      timeStamps = as.Date(gridTime$Var2), extremeValues = Y,
      pdf = pdf
    )
  }
  # ggplot is really too slow


  # rgb.palette=colorRampPalette(c("F6FF33","#ffeda0", "#fed976", "#feb24c","#fd8d3c","#fc4e2a",
  #                                "#e31a1c","#bd0026","#800026",'#2C110B'),interpolate="linear",bias=1)

  rgb.palette <- colorRampPalette(rev(c(
    "#2C110B", "#8B0000", "#FF0000", "#FF4500", "#FFA500",
    "#FFD700", "#FFFF00", "#FFFFE0"
  )), interpolate = "linear", bias = 1)

  # datap <- data.frame(timeStamps = gridTime$Var2, extremeValues = Y,
  #                     pdf = pdf)

  # time breaks
  tbi <- round(lubridate::year(minTS) / 5) * 5
  tbf <- round(lubridate::year(maxTS) / 5) * 5

  ttbi <- as.Date(paste(tbi, "-01-01", sep = ""))
  ttbf <- as.Date(paste(tbf, "-01-01", sep = ""))

  ylims <- c(max(min(serix, na.rm = T) - 1, 0), min(max(Y)))
  bY <- round((ylims[2] - ylims[1]) / 100) * 10
  ey <- ylims[1]
  if (trans == "inv") {
    peakplot$value <- 1 / peakplot$value
    ylims <- c(max(min(1 / serix, na.rm = T) - 1, 0), max(Ybis))
    serie$serio <- 1 / serie$serix
    ylims <- c(max(min(serie$serio, na.rm = T) - 1, 0), min(max(1 / Y)))
    ey <- ylims[2]
  } else if (trans == "rev") {
    peakplot$value <- -peakplot$value
    serie$serio <- -serie$serix
    ylims <- c(max(min(serie$serio, na.rm = T) - sd(peakplot$value, na.rm = T), 0), min(max(-Y)))
    print(ylims)
    ey <- ylims[2]
  } else if (trans == "lninv") {
    peakplot$value <- 1 / exp(peakplot$value)
    serie$serio <- 1 / exp(serie$serix)
    ylims <- c(max(min(serie$serio, na.rm = T) - 1, 0), min(max(1 / exp(Y))))
    ey <- ylims[2]
  } else {
    serie$serio <- serie$serix
  }

  datap$npdf <- datap$pdf / max(datap$pdf)
  plo <- ggplot(datap, aes(x = as.Date(timeStamps), y = extremeValues)) +
    geom_segment(data = serie, aes(x = timeStamps, xend = timeStamps, y = ey, yend = serio), col = "black", alpha = .5) +
    geom_raster(alpha = (datap$npdf)^0.2, aes(fill = pdf), interpolate = TRUE) +
    scale_fill_gradientn(
      colours = rgb.palette(100),
      n.breaks = 10, guide = "coloursteps", trans = scales::modulus_trans(0.7), na.value = "transparent"
    ) +
    ggtitle(paste0("GPD - ", args$Title)) +
    geom_point(
      data = peakplot, aes(x = as.Date(time), y = value), shape = 21, fill = "white",
      color = "black", size = 2, stroke = 2
    ) +
    scale_y_continuous(
      n.breaks = 10,
      limit = c(ylims[1], ylims[2]), args$ylabel, expand = c(0, 0)
    ) +
    scale_x_date(
      labels = scales::date_format("%Y"), args$xlabel, breaks = seq(ttbi, ttbf, by = paste0(tic, " years")),
      minor_breaks = scales::date_breaks(paste0(tac, " months")), expand = c(0, 0), limits = c(minTS, maxTS)
    ) +
    annotate("label", x = maxTS - 36 * xy, y = 0.95 * ylims[2], label = paste0("\U03B5 = ", as.character(round(epsilon, 3))), size = 8) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = args$axisFontSize - 6),
      panel.grid.major.x = element_line(color = "black"),
      axis.text.y = element_text(size = args$axisFontSize),
      axis.title.x = element_text(size = args$labelFontSize),
      axis.title.y = element_text(size = args$labelFontSize)
    ) +
    labs(x = args$xlabel, y = args$ylabel)
  return(plo)
}

#' Generate a GEV (Generalized Extreme Value) plot with raster image
#'
#' This function generates a GEV plot with a raster image using the Generalized Extreme Value (GEV) distribution.
#'
#' @param Y A vector of extreme values.
#' @param timeStamps A vector of timestamps corresponding to the extreme values.
#' @param serix The y-value at which to draw a horizontal line on the plot.
#' @param epsilon A numeric value representing the shape parameter of the GEV distribution.
#' @param sigma A vector of scale parameters corresponding to the timestamps.
#' @param mu A vector of location parameters corresponding to the timestamps.
#' @param returnPeriodInDts The return period in decimal time steps.
#' @param maxObs A data frame containing the maximum observations.
#' @param mode A character string indicating the mode of the plot. Possible values are "inv" (inverse) and "normal".
#' @param ... Additional arguments to customize the plot.
#'
#' @return A ggplot object representing the GEV plot with a raster image.
#'
#' @examples
#' # Generate a GEV plot with raster image
#' tsEvaPlotGEVImageSc(
#'   Y = c(1, 2, 3), timeStamps = c("2021-01-01", "2021-01-02", "2021-01-03"),
#'   serix = 2, epsilon = 0.1, sigma = c(0.5, 0.6, 0.7),
#'   mu = c(1, 2, 3), returnPeriodInDts = 1, maxObs = data.frame(time = "2021-01-01", value = 2),
#'   mode = "normal"
#' )
#'
#' @import ggplot2 scales
#' @importFrom texmex dgev
#'
#' @export
tsEvaPlotGEVImageSc <- function(Y, timeStamps, serix, epsilon, sigma, mu, returnPeriodInDts, maxObs, mode, varargin) {
  avgYearLength <- 365.2425
  nyears <- as.numeric(round((max(timeStamps) - min(timeStamps)) / avgYearLength))
  nelmPerYear <- length(timeStamps) / nyears

  args <- list(
    nPlottedTimesByYear = min(5, round(nelmPerYear)),
    ylabel = "levels (mm)",
    xlabel = "Date",
    minYear = 1950,
    maxYear = 2021,
    axisFontSize = 20,
    labelFontSize = 22,
    Title = ""
  )
  args <- tsEasyParseNamedArgs(varargin, args)

  minTS <- as.Date(paste(args$minYear, "-01-01", sep = ""))
  maxTS <- as.Date(paste(args$maxYear, "-01-01", sep = ""))
  sigma <- sigma[(timeStamps >= minTS) & (timeStamps <= maxTS)]
  mu <- mu[(timeStamps >= minTS) & (timeStamps <= maxTS)]
  timeStamps <- timeStamps[(timeStamps >= minTS) & (timeStamps <= maxTS)]


  L <- length(timeStamps)
  minTS <- timeStamps[1]
  maxTS <- timeStamps[length(timeStamps)]
  npdf <- as.numeric(ceiling(((maxTS - minTS) / avgYearLength) * args$nPlottedTimesByYear))
  navg <- floor(L / npdf)
  navg <- 1

  plotSLength <- npdf * navg
  timeStamps_plot <- seq(minTS, maxTS, length.out = plotSLength)

  if (length(epsilon) == 1) {
    epsilon0 <- rep(epsilon, npdf)
  } else {
    epsilon_ <- rep(NA, npdf * navg)
    epsilon_[1:L] <- epsilon
    epsilonMtx <- matrix(epsilon_, navg, ncol = npdf)
    epsilon0 <- apply(epsilonMtx, 2, mean, na.rm = TRUE)
  }

  sigma_ <- approx(timeStamps, sigma, timeStamps_plot)$y
  sigmaMtx <- matrix(sigma_, navg, ncol = npdf)
  if (length(sigmaMtx) > 1) {
    sigma0 <- apply(sigmaMtx, 2, mean, na.rm = TRUE)
  } else {
    sigma0 <- sigmaMtx
  }

  mu_ <- approx(timeStamps, mu, timeStamps_plot)$y
  muMtx <- matrix(mu_, navg, ncol = npdf)
  if (length(muMtx) > 1) {
    mu0 <- apply(muMtx, 2, mean, na.rm = TRUE)
  } else {
    mu0 <- muMtx
  }

  epsilon0 <- as.numeric(epsilon0)
  sigma0 <- as.numeric(sigma0)
  mu0 <- as.numeric(mu0)
  Ys <- Y

  gridEps <- expand.grid(Y, epsilon0)
  gridSig <- expand.grid(Y, sigma0)
  gridMu <- expand.grid(Y, mu0)
  gridTime <- expand.grid(Y, timeStamps_plot)

  pdf <- texmex::dgev(Y, gridMu$Var2, gridSig$Var2, gridEps$Var2)

  if (mode == "inv") {
    Ybis <- seq(min(1 / Y), max(1 / Y), length.out = length(Y))
    pdf <- texmex::dgev(1 / Ybis, gridMu$Var2, gridSig$Var2, gridEps$Var2)
    datap <- data.frame(
      timeStamps = gridTime$Var2, extremeValues = Ybis,
      pdf = pdf
    )
  } else {
    pdf <- texmex::dgev(Y, gridMu$Var2, gridSig$Var2, gridEps$Var2)
    datap <- data.frame(
      timeStamps = gridTime$Var2, extremeValues = Y,
      pdf = pdf
    )
  }
  # ggplot is really too slow


  # rgb.palette=colorRampPalette(c("#EBFE00","#F6FF33","#ffeda0", "#fed976", "#feb24c","#fd8d3c","#fc4e2a",
  #                                "#e31a1c","#bd0026","#800026",'#2C110B'),interpolate="linear",bias=1)

  rgb.palette <- colorRampPalette(rev(c(
    "#2C110B", "#8B0000", "#FF0000", "#FFA500",
    "#FFF700", "#FFFFF6"
  )), interpolate = "linear", bias = 1)


  # time breaks
  tbi <- round(lubridate::year(minTS) / 5) * 5
  tbf <- round(lubridate::year(maxTS) / 5) * 5

  ttbi <- as.Date(paste(tbi, "-01-01", sep = ""))
  ttbf <- as.Date(paste(tbf, "-01-01", sep = ""))
  xy <- (tbf - tbi)
  if (xy >= 20) {
    tic <- 5
    tac <- 12
  }
  if (xy < 10) {
    tic <- 1
    tac <- 1
    tbi <- lubridate::year(minTS)
  }
  if (xy >= 10 & xy < 20) {
    tic <- 2
    tac <- 6
  }

  ttbi <- as.Date(paste(tbi, "-01-01", sep = ""))
  ttbf <- as.Date(paste(tbf, "-01-01", sep = ""))
  # correction to get a more usefull graph

  ylims <- c(max(min(Y), 0), min(max(Ys)))



  if (mode == "inv") {
    maxObs$value <- 1 / maxObs$value
    ylims <- c(max(min(Ybis) - 0.5, 0), max(Ybis) + 0.5)
  }

  datapsub <- datap
  datapsub$npdf <- datapsub$pdf / max(datapsub$pdf, na.rm = T)

  plo <- ggplot(datapsub, aes(x = timeStamps, y = extremeValues)) +
    geom_raster(alpha = (datapsub$npdf)^0.2, aes(fill = pdf), interpolate = T) +
    scale_fill_gradientn(
      colours = rgb.palette(100),
      n.breaks = 10, guide = "coloursteps", trans = scales::modulus_trans(1.1), na.value = "transparent"
    ) +
    ggtitle(paste0("GEV - ", args$Title)) +
    geom_point(
      data = maxObs, aes(x = as.Date(time), y = value), shape = 21, fill = "white",
      color = "black", size = 1, stroke = 2
    ) +
    scale_y_continuous(
      n.breaks = 10,
      limit = c(ylims[1], ylims[2]), args$ylabel, expand = c(0, 0)
    ) +
    annotate("label", x = maxTS - 36 * xy, y = 0.9 * ylims[2], label = paste0("\U03B5 = ", as.character(round(epsilon, 3))), size = 8) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = args$axisFontSize),
      panel.grid.major.x = element_line(color = "black"),
      axis.text.y = element_text(size = args$axisFontSize),
      axis.title.x = element_text(size = args$labelFontSize),
      axis.title.y = element_text(size = args$labelFontSize)
    ) +
    labs(x = args$xlabel, y = args$ylabel)
  plo
  return(plo)
}




#' Convert Transformed Time Series to Stationary Series and Plot
#'
#' This function takes the parameters of a non-stationary time series evaluation,
#'  along with the transformed stationary data, and plots the converted stationary series.
#'
#' @param nonStationaryEvaParams A list of parameters for non-stationary time series evaluation.
#' @param stationaryTransformData A list containing the transformed stationary data.
#' @param ... Additional arguments to be passed to the \code{\link{tsEvaPlotTransfToStat}} function.
#'
#' @import ggplot2
#' @seealso \code{\link{tsEvaPlotTransfToStat}}
#' @return The plot object representing the converted stationary series.
#'
#' @examples
#' # Example usage
#' nonStationaryParams <- list(param1 = 1, param2 = 2)
#' stationaryData <- list(timeStamps = c(1, 2, 3), stationarySeries = c(0.1, 0.2, 0.3), statSer3Mom = c(0.01, 0.02, 0.03), statSer4Mom = c(0.001, 0.002, 0.003))
#' tsEvaPlotTransfToStatFromAnalysisObj(nonStationaryParams, stationaryData)
#'
#' @export
tsEvaPlotTransfToStatFromAnalysisObj <- function(nonStationaryEvaParams, stationaryTransformData, ...) {
  varargin <- NULL
  varargin <- list(...)
  # print(varargin)
  # Extract timeStamps, series, mean, stdDev, 3rd moment, 4th moment
  timeStamps <- stationaryTransformData$timeStamps
  statSeries <- stationaryTransformData$stationarySeries
  srsmean <- rep(0, length(statSeries))
  stdDev <- rep(1, length(statSeries))
  st3mom <- stationaryTransformData$statSer3Mom
  st4mom <- stationaryTransformData$statSer4Mom

  plotbg <- tsEvaPlotTransfToStat(timeStamps, statSeries, srsmean, stdDev, st3mom, st4mom, varargin)

  print(plotbg)

  return(plotbg)
}


#' Plot Time Series with Trend and Standard Deviation
#'
#' This function plots a time series along with its trend and standard deviation.
#'
#' @param nonStationaryEvaParams The non-stationary evaluation parameters.
#' @param stationaryTransformData The stationary transformed data.
#' @param mode The mode of the plot (optional).
#' @param ... Additional arguments to customize the plot (optional).
#'
#' @import ggplot2
#' @seealso \code{\link{trans}}
#' @return A ggplot object representing the plot.
#'
#' @examples
#' # Example usage
#' tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStationaryEvaParams, stationaryTransformData, mode = "inv")
#' @export
tsEvaPlotSeriesTrendStdDevFromAnalyisObj <- function(nonStationaryEvaParams, stationaryTransformData, mode, ...) {
  varargin <- list(...)
  args <- list(
    plotPercentile = -1,
    ylabel = "Values",
    xlabel = "Date",
    minYear = 1950,
    maxYear = 2020,
    axisFontSize = 20,
    labelFontSize = 22
  )
  args <- tsEasyParseNamedArgs(varargin, args)

  minTS <- as.Date(paste(args$minYear, "-01-01", sep = ""))
  maxTS <- as.Date(paste(args$maxYear, "-01-01", sep = ""))
  args <- tsEasyParseNamedArgs(varargin, args)
  plotPercentile <- args$plotPercentile

  timeStamps <- stationaryTransformData$timeStamps
  series <- stationaryTransformData$nonStatSeries
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  if (dt == 1) {
    trend <- stationaryTransformData$trendSeries
    stdDev <- stationaryTransformData$stdDevSeries
  } else {
    trend <- stationaryTransformData$trendSeriesOr
    stdDev <- stationaryTransformData$stdDevSeriesOr
  }
  if ("statsTimeStamps" %in% names(stationaryTransformData)) {
    statsTimeStamps <- stationaryTransformData$statsTimeStamps
  } else {
    statsTimeStamps <- timeStamps
  }

  # Create a data frame with the time stamps, series, trend, and stdDev
  data <- data.frame(timeStamps = timeStamps, series = series, trend = trend, stdDev = stdDev)
  data$infCI <- data$trend - data$stdDev
  data$supCI <- data$trend + data$stdDev
  # time breaks
  tbi <- round(lubridate::year(minTS) / 5) * 5
  tbf <- round(lubridate::year(maxTS) / 5) * 5

  ttbi <- as.Date(paste(tbi, "-01-01", sep = ""))
  ttbf <- as.Date(paste(tbf, "-01-01", sep = ""))
  if (mode == "inv") {
    data$series <- 1 / data$series
    data$infCI <- 1 / data$infCI
    data$supCI <- 1 / data$supCI
    data$trend <- 1 / data$trend
  }
  # Create the plot
  plotz <- ggplot(data, aes(x = as.Date(timeStamps), y = series)) +
    geom_line(aes(color = "Series"), size = 1) +
    geom_ribbon(aes(ymin = infCI, ymax = supCI), fill = "lightgreen", alpha = 0.6) +
    geom_line(aes(y = infCI, color = "Std.Dev"), size = 1, lty = 2) +
    geom_line(aes(y = supCI, color = "Std.Dev"), size = 1, lty = 2) +
    geom_line(aes(y = trend, color = "Trend"), size = 1) +
    ggtitle("Trend") +
    scale_color_manual(name = "", values = c("Trend" = "black", "Std.Dev" = "darkgreen", "Series" = "red")) +
    scale_y_continuous(
      n.breaks = 5, args$ylabel
    ) +
    scale_x_date(
      labels = scales::date_format("%Y"), args$xlabel, breaks = seq(ttbi, ttbf, by = "5 years"), expand = c(0, 0)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = args$axisFontSize),
      axis.text.y = element_text(size = args$axisFontSize),
      axis.title.x = element_text(size = args$labelFontSize),
      axis.title.y = element_text(size = args$labelFontSize),
      legend.justification = c(1.1, 1.1), legend.position = c(1, 1),
      legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
      legend.title = element_blank()
    )
  plotz
  return(plotz)
}


#' tsEvaPlotTransfToStat Function
#'
#' This function creates a line plot of time series data along with statistical measures.
#'
#' @param timeStamps A vector of time stamps for the data points.
#' @param statSeries A vector of the main time series data.
#' @param srsmean A vector of the mean values for each time stamp.
#' @param stdDev A vector of the standard deviation values for each time stamp.
#' @param st3mom A vector of the third moment values for each time stamp.
#' @param st4mom A vector of the fourth moment values for each time stamp.
#' @param varargin Additional optional arguments to customize the plot.
#'
#' @import ggplot2
#' @return A ggplot object representing the line plot.
#'
#' @examples
#' # Example usage of tsEvaPlotTransfToStat function
#' data <- data.frame(
#'   timeStamps = c("2021-01-01", "2021-01-02", "2021-01-03"),
#'   statSeries = c(1, 2, 3),
#'   srsmean = c(1.5, 2.5, 3.5),
#'   stdDev = c(0.5, 0.6, 0.7),
#'   st3mom = c(0.1, 0.2, 0.3),
#'   st4mom = c(0.01, 0.02, 0.03)
#' )
#' plot <- tsEvaPlotTransfToStat(data$timeStamps, data$statSeries, data$srsmean, data$stdDev, data$st3mom, data$st4mom)
#' print(plot)
#'
#' @export
tsEvaPlotTransfToStat <- function(timeStamps, statSeries, srsmean, stdDev, st3mom, st4mom, varargin) {
  data <- data.frame(timeStamps, statSeries, srsmean, stdDev, st3mom, st4mom)
  names(data) <- c("timeStamps", "statSeries", "srsmean", "stdDev", "thirdMom", "fourthMom")


  avgYearLength <- 365.2425
  nyears <- as.numeric(round((max(timeStamps) - min(timeStamps)) / avgYearLength))
  nelmPerYear <- length(timeStamps) / nyears

  args <- list(
    nPlottedTimesByYear = min(360, round(nelmPerYear)),
    ylabel = " ",
    xlabel = "Date",
    minYear = 1979,
    maxYear = 2019,
    axisFontSize = 20,
    labelFontSize = 22,
    legendFontSize = 18
  )

  args <- tsEasyParseNamedArgs(varargin, args)


  # time breaks
  minTS <- data$timeStamps[1]
  maxTS <- data$timeStamps[length(timeStamps)]

  tbi <- round(lubridate::year(minTS) / 5) * 5
  tbf <- round(lubridate::year(maxTS) / 5) * 5

  ttbi <- as.Date(paste(tbi, "-01-01", sep = ""))
  ttbf <- as.Date(paste(tbf, "-01-01", sep = ""))

  ylims <- c(min(statSeries, na.rm = T), max(max(statSeries, na.rm = T), max(st3mom), max(st4mom)))

  elplot <- ggplot(data, aes(x = timeStamps)) +
    geom_line(aes(y = statSeries, color = "Normalized series"), color = "royalblue") +
    geom_line(aes(y = srsmean, color = "Mean"), linetype = "dashed", size = 1.5) +
    geom_line(aes(y = stdDev, color = "Std.Dev"), linetype = "dashed", size = 1.5) +
    geom_line(aes(y = thirdMom, color = "Skewness"), color = "purple", size = 1.5) +
    geom_line(aes(y = fourthMom, color = "Kurtosis"), color = "darkgreen", size = 1.5) +
    scale_color_manual(name = "", values = c("Normalized series" = "royalblue", "Mean" = "black", "Std.Dev" = "red", "Skewness" = "purple", "Kurtosis" = "darkgreen")) +
    scale_y_continuous(
      breaks = seq(round(ylims[1] / 10) * 10, ylims[2], by = 5),
      limit = c(ylims[1], ylims[2] + 1), args$ylabel
    ) +
    scale_x_date(
      labels = scales::date_format("%Y"), args$xlabel, breaks = seq(ttbi, ttbf, by = "5 years"), expand = c(0, 0)
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = args$axisFontSize),
      axis.title = element_text(size = args$axisFontSize),
      legend.text = element_text(size = args$legendFontSize),
      legend.justification = c(1.1, 1.1), legend.position = c(1, 1),
      legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
      legend.title = element_blank()
    )

  return(elplot)
}
