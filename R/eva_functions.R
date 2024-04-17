# EVA functions to compute relturn levels / returtn periods from Nonstationary EVD parameters-----------------

#' tsEvaComputeReturnPeriodsGEV
#'
#' This function computes the return levels and return periods for a Generalized Extreme Value (GEV) distribution,
#' given the GEV parameters and their standard error. The return levels represent the values of annual maxima
#' with a certain probability, while the return periods indicate the average time between
#' exceedances of those threshold values.
#'
#' @param epsilon The shape parameter of the GEV distribution.
#' @param sigma The scale parameter of the GEV distribution.
#' @param mu The location parameter of the GEV distribution.
#' @param BlockMax A vector containing the block maxima data.
#' @param MaxID An identifier for the maximum value used in the computation.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{GevPseudo}: A matrix of pseudo observations obtained from the GEV distribution for each annual extreme at every time step.
#' \item \code{returnPeriods}: A matrix of return periods corresponding to the pseudo observations.
#' \item \code{PseudoObs}: The pseudo observation corresponding to the maximum value used in the computation.
#' }
#'
#' @seealso \code{\link{empdis}}
#'
#' @family MyPackage functions
#'
#' @docType methods
#' @name tsEvaComputeReturnPeriodsGEV
#' @export
tsEvaComputeReturnPeriodsGEV <- function(epsilon, sigma, mu, BlockMax, MaxID) {
  # tsEvaComputeReturnLevelsGEV:
  # returns the return levels given the gev parameters and their standard
  # error.
  # The parameter returnPeriodsInDts contains the return period expressed in
  # a time unit that corresponds to the size of the time segments where we
  # are evaluating the maxima. For example, if we are working on yearly
  # maxima, returnPeriodsInDts must be expressed in years. If we are working
  # on monthly maxima returnPeriodsInDts must be expressed in months.


  nyr <- length(BlockMax)
  # Compute the empirical return period:
  empval <- empdis(BlockMax, nyr)
  my <- match(BlockMax, empval$Q)
  PseudoObs <- empval[my, ]

  # uniforming dimensions of yp, sigma, mu
  npars <- length(sigma)
  nt <- length(BlockMax)
  Bmax <- matrix(BlockMax, nrow = npars, ncol = nt, byrow = TRUE)
  sigma_ <- matrix(sigma, nrow = npars, ncol = nt, byrow = FALSE)
  mu_ <- matrix(mu, nrow = npars, ncol = nt, byrow = FALSE)

  # estimating the return levels
  qx <- 1 - exp(-(1 + epsilon * (Bmax - mu_) / sigma_)^(-1 / epsilon))
  GevPseudo <- qx
  returnPeriods <- 1 / qx

  return(list(GevPseudo = GevPseudo, returnPeriods = returnPeriods, PseudoObs = PseudoObs))
}

#' tsEvaComputeReturnPeriodsGPD
#'
#' This function computes the return periods and pseudo observations for a
#' Generalized Pareto Distribution (GPD), given the GPD parameters,
#' threshold, peaks data, and sample time horizon.
#'
#' @param epsilon The shape parameter of the GPD.
#' @param sigma The scale parameter of the GPD.
#' @param threshold The threshold value for the GPD.
#' @param peaks A vector containing the peak values.
#' @param nPeaks The number of peak values.
#' @param peaksID An identifier for each peak value.
#' @param sampleTimeHorizon The time horizon of the sample.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{GpdPseudo}: A matrix of pseudo observations obtained from the GPD for each peak value at every time step.
#'     \item \code{returnPeriods}: A matrix of return periods corresponding to the pseudo observations.
#'     \item \code{PseudoObs}: A data frame containing the pseudo observations and their corresponding identifiers.
#'   }
#'
#' @references
#' Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L., Besio, G., and Alfieri, L. (2016). The transformed-stationary approach: a generic and simplified methodology for non-stationary extreme value analysis. \emph{Hydrology and Earth System Sciences}, \strong{20}(9), 3527-3547. doi:10.5194/hess-20-3527-2016.
#'
#' @keywords function, computation, GPD, return periods, pseudo observations
#' @seealso \code{\link{empdis}}
#'
#'
#' @docType methods
#' @name tsEvaComputeReturnPeriodsGPD
#' @export
tsEvaComputeReturnPeriodsGPD <- function(epsilon, sigma, threshold, peaks, nPeaks, peaksID, sampleTimeHorizon) {
  X0 <- nPeaks / sampleTimeHorizon
  nyr <- sampleTimeHorizon
  peaksj <- jitter(peaks, 2)
  peakAndId <- data.frame(peaksID, peaksj)
  peaksjid <- peaksID[order(peaksj)]

  # Compute the empirical return period:
  empval <- empdis(peaksj, nyr)
  PseudoObs <- empval
  PseudoObs$ID <- peaksjid
  PseudoObs <- PseudoObs[order(PseudoObs$ID), ]

  # Uniform dimensions of variables
  npars <- length(sigma)
  nt <- length(peaks)
  Peakx <- matrix(PseudoObs$Q, nrow = npars, ncol = nt, byrow = TRUE)
  sigma_ <- matrix(sigma, nrow = npars, ncol = nt, byrow = FALSE)
  threshold_ <- matrix(threshold, nrow = npars, ncol = nt, byrow = FALSE)

  # Estimate the return levels
  # Here I have all the pseudo observations of every annual extreme for every timestep
  qx <- (((1 + epsilon * (Peakx - threshold_) / sigma_)^(-1 / epsilon)))
  GpdPseudo <- qx
  GpdPseudo[which(GpdPseudo < 0)] <- 0
  returnPeriods <- 1 / (X0 * qx)

  return(list(GpdPseudo = GpdPseudo, returnPeriods = returnPeriods, PseudoObs = PseudoObs))
}

#' tsEvaComputeReturnLevelsGEV
#'
#' This function calculates the return levels for a Generalized Extreme Value (GEV)
#' distribution given the GEV parameters and their standard errors.
#' The return periods are specified in a time unit that corresponds
#' to the size of the time segments for evaluating the maxima.
#'
#' @param epsilon The shape parameter of the GEV distribution.
#' @param sigma The scale parameter of the GEV distribution.
#' @param mu The location parameter of the GEV distribution.
#' @param epsilonStdErr The standard error of the shape parameter.
#' @param sigmaStdErr The standard error of the scale parameter.
#' @param muStdErr The standard error of the location parameter.
#' @param returnPeriodsInDts The return periods expressed in the corresponding time unit (e.g., years for yearly maxima).
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{returnLevels}: A matrix of return levels corresponding to the specified return periods.
#'     \item \code{returnLevelsErr}: A matrix of standard errors for the return levels.
#'   }
#'
#' @references
#' Stuart Coles (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}. Springer.
#' Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L.,
#' Besio, G., and Alfieri, L. (2016). The transformed-stationary approach:
#' a generic and simplified methodology for non-stationary extreme value analysis.
#' \emph{Hydrology and Earth System Sciences}, \strong{20}(9), 3527-3547.
#'  doi:10.5194/hess-20-3527-2016.
#'
#' @keywords function, computation, GEV, return levels, standard errors
#' @seealso \code{\link{empdis}}
#'
#'
#' @docType methods
#' @name tsEvaComputeReturnLevelsGEV
#' @export
tsEvaComputeReturnLevelsGEV <- function(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriodsInDts) {
  # tsEvaComputeReturnLevelsGEV:
  # returns the return levels given the gev parameters and their standard
  # error.
  # The parameter returnPeriodsInDts contains the return period expressed in
  # a time unit that corresponds to the size of the time segments where we
  # are evaluating the maxima. For example, if we are working on yearly
  # maxima, returnPeriodsInDts must be expressed in years. If we are working
  # on monthly maxima returnPeriodsInDts must be expressed in months.

  returnPeriodsInDts <- as.vector(returnPeriodsInDts)
  returnPeriodsInDtsSize <- length(returnPeriodsInDts)
  if (returnPeriodsInDtsSize[1] > 1) {
    returnPeriodsInDts <- as.vector(t(returnPeriodsInDts))
  }
  # reference: Stuart Coles 2001, pag 49.
  yp <- -log(1 - 1 / returnPeriodsInDts)

  # uniforming dimensions of yp, sigma, mu
  npars <- length(sigma)
  nt <- length(returnPeriodsInDts)
  yp <- matrix(yp, nrow = npars, ncol = nt, byrow = TRUE)
  sigma_ <- matrix(sigma, nrow = npars, ncol = nt, byrow = TRUE)
  sigmaStdErr_ <- matrix(sigmaStdErr, nrow = npars, ncol = nt, byrow = TRUE)
  mu_ <- matrix(mu, nrow = npars, ncol = nt, byrow = TRUE)
  muStdErr_ <- matrix(muStdErr, nrow = npars, ncol = nt, byrow = TRUE)
  if (epsilon != 0) {
    # estimating the return levels
    returnLevels <- mu_ - sigma_ / epsilon * (1 - yp^(-epsilon))
    ## estimating the error
    # estimating the differential of returnLevels to the parameters
    dxm_mu <- 1
    dxm_sigma <- 1 / epsilon * (1 - yp^(-epsilon))
    dxm_epsilon <- sigma_ / (epsilon^2) * (1 - (yp)^(-epsilon)) - sigma_ / epsilon * log(yp) * yp^(-epsilon)
    returnLevelsErr <- sqrt((dxm_mu * muStdErr_)^2 + (dxm_sigma * sigmaStdErr_)^2 + (dxm_epsilon * epsilonStdErr)^2)
  } else {
    returnLevels <- mu_ - sigma_ * log(yp)
    ## estimating the error
    # estimating the differential of returnLevels to the parameter
    dxm_u <- rep(1, length(mu_))
    dxm_sigma <- log(yp)
    returnLevelsErr <- sqrt((dxm_u * muStdErr_)^2 + (dxm_sigma * sigmaStdErr_)^2)
  }
  return(list(returnLevels = returnLevels, returnLevelsErr = returnLevelsErr))
}

#' tsEvaComputeReturnLevelsGEVFromAnalysisObj
#'
#' This function calculates the return levels for a Generalized Extreme Value (GEV) distribution using the parameters obtained from a non-stationary extreme value analysis. It supports non-stationary analysis by considering different parameters for each time index.
#'
#' @param nonStationaryEvaParams The parameters obtained from a non-stationary extreme value analysis.
#' @param returnPeriodsInYears The return periods expressed in years.
#' @param ... Additional arguments.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{returnLevels}: A matrix of return levels corresponding to the specified return periods.
#'     \item \code{returnLevelsErr}: A matrix of standard errors for the return levels.
#'     \item \code{returnLevelsErrFit}: A matrix of standard errors for the return levels obtained from fitting the non-stationary model.
#'     \item \code{returnLevelsErrTransf}: A matrix of standard errors for the return levels obtained from the transformed data.
#'   }
#'
#' @seealso \code{\link{tsEvaComputeReturnLevelsGEV}}
#'
#' @family MyPackage functions
#'
#' @docType methods
#' @name tsEvaComputeReturnLevelsGEVFromAnalysisObj
#' @export
tsEvaComputeReturnLevelsGEVFromAnalysisObj <- function(nonStationaryEvaParams, returnPeriodsInYears, ...) {
  args <- list(timeIndex = -1)
  varargin <- list(...)
  args <- tsEasyParseNamedArgs(varargin, args)
  timeIndex <- args$timeIndex

  epsilon <- nonStationaryEvaParams[[1]]$parameters$epsilon
  epsilonStdErr <- nonStationaryEvaParams[[1]]$paramErr$epsilonErr
  epsilonStdErrFit <- epsilonStdErr
  epsilonStdErrTransf <- 0
  marker <- nonStationaryEvaParams[[1]]$paramErr$sigmaErrTransf
  nonStationary <- exists("marker")
  if (timeIndex > 0) {
    sigma <- nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex]
    mu <- nonStationaryEvaParams[[1]]$parameters$mu[timeIndex]
    sigmaStdErr <- nonStationaryEvaParams[[1]]$paramErr$sigmaErr[timeIndex]
    sigmaStdErrFit <- nonStationaryEvaParams[[1]]$paramErr$sigmaErrFit[timeIndex]
    sigmaStdErrTransf <- nonStationaryEvaParams[[1]]$paramErr$sigmaErrTransf[timeIndex]
    muStdErr <- nonStationaryEvaParams$paramErr[[1]]$muErr[timeIndex]
    muStdErrFit <- nonStationaryEvaParams[[1]]$paramErr$muErrFit[timeIndex]
    muStdErrTransf <- nonStationaryEvaParams[[1]]$paramErr$muErrTransf[timeIndex]
  } else {
    sigma <- nonStationaryEvaParams[[1]]$parameters$sigma
    mu <- nonStationaryEvaParams[[1]]$parameters$mu
    sigmaStdErr <- nonStationaryEvaParams[[1]]$paramErr$sigmaErr
    muStdErr <- nonStationaryEvaParams[[1]]$paramErr$muErr
    if (nonStationary) {
      sigmaStdErrFit <- nonStationaryEvaParams[[1]]$paramErr$sigmaErrFit
      sigmaStdErrTransf <- nonStationaryEvaParams[[1]]$paramErr$sigmaErrTransf
      muStdErrFit <- nonStationaryEvaParams[[1]]$paramErr$muErrFit
      muStdErrTransf <- nonStationaryEvaParams[[1]]$paramErr$muErrTransf
    }
  }

  returnLevels <- tsEvaComputeReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriodsInYears)
  returnLevelsErr <- returnLevels$returnLevelsErr
  if (nonStationary) {
    returnLevelsErrFit <- tsEvaComputeReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErrFit, sigmaStdErrFit, muStdErrFit, returnPeriodsInYears)$returnLevelsErr
    returnLevelsErrTransf <- tsEvaComputeReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErrTransf, sigmaStdErrTransf, muStdErrTransf, returnPeriodsInYears)$returnLevelsErr
  } else {
    returnLevelsErrFit <- returnLevelsErr
    returnLevelsErrTransf <- rep(0, length(returnLevelsErr))
  }

  return(list(returnLevels = returnLevels$returnLevels, returnLevelsErr = returnLevelsErr, returnLevelsErrFit = returnLevelsErrFit, returnLevelsErrTransf = returnLevelsErrTransf))
}

#' tsEvaComputeReturnLevelsGPD
#'
#' This function calculates the return levels for a Generalized Pareto Distribution (GPD) using the parameters of the distribution and their standard errors.
#'
#' @param epsilon The shape parameter of the GPD.
#' @param sigma The scale parameter of the GPD.
#' @param threshold The threshold parameter of the GPD.
#' @param epsilonStdErr The standard error of the shape parameter.
#' @param sigmaStdErr The standard error of the scale parameter.
#' @param thresholdStdErr The standard error of the threshold parameter.
#' @param nPeaks The number of peaks used to estimate the parameters.
#' @param sampleTimeHorizon The time horizon of the sample in the same units as the return periods (e.g., years).
#' @param returnPeriods The return periods for which to compute the return levels.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{returnLevels}: A vector of return levels corresponding to the specified return periods.
#'     \item \code{returnLevelsErr}: A vector of standard errors for the return levels.
#'   }
#' @references
#' Stuart Coles (2001). \emph{An Introduction to Statistical Modeling of Extreme Values}. Springer.
#' @details
#' sampleTimeHorizon and returnPeriods must be in the same units, e.g. years
#'
#' @seealso \code{\link{tsEasyParseNamedArgs}}
#'
#' @family MyPackage functions
#'
#' @docType methods
#' @name tsEvaComputeReturnLevelsGPD
#' @export
tsEvaComputeReturnLevelsGPD <- function(epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, sampleTimeHorizon, returnPeriods) {

  X0 <- nPeaks / sampleTimeHorizon
  XX <- X0 * returnPeriods
  npars <- length(sigma)
  nt <- length(returnPeriods)
  XX_ <- rep(XX, each = npars)
  sigma_ <- rep(sigma, times = nt)
  sigmaStdErr_ <- rep(sigmaStdErr, times = nt)
  threshold_ <- rep(threshold, times = nt)
  thresholdStdErr_ <- rep(thresholdStdErr, times = nt)

  if (epsilon != 0) {
    # estimating the return levels
    returnLevels <- threshold_ + sigma_ / epsilon * ((XX_)^epsilon - 1)
    # estimating the error
    # estimating the differential of returnLevels to the parameters
    # !! ASSUMING NON ZERO ERROR ON THE THRESHOLD AND 0 ERROR ON THE PERCENTILE.
    # THIS IS NOT COMPLETELY CORRECT BECAUSE THE PERCENTILE DEPENDS ON THE
    # THRESHOLD AND HAS THEREFORE AN ERROR RELATED TO THAT OF THE
    # THRESHOLD
    dxm_u <- 1
    dxm_sigma <- 1 / epsilon * (XX_^epsilon - 1)
    dxm_epsilon <- -sigma_ / epsilon^2 * ((XX_)^epsilon - 1) + sigma_ / epsilon * log(XX_) * XX_^epsilon

    returnLevelsErr <- sqrt((dxm_u * thresholdStdErr_)^2 + (dxm_sigma * sigmaStdErr_)^2 + (dxm_epsilon * epsilonStdErr)^2)
    #
  } else {
    returnLevels <- threshold_ + sigma_ * log(XX_)
    dxm_u <- 1
    dxm_sigma <- log(XX_)

    returnLevelsErr <- ((dxm_u * thresholdStdErr_)^2 + (dxm_sigma * sigmaStdErr_)^2)^.5
  }
  return(list(returnLevels = returnLevels, returnLevelsErr = returnLevelsErr))
}

#' tsEvaComputeReturnLevelsGPDFromAnalysisObj
#'
#' \code{tsEvaComputeReturnLevelsGPDFromAnalysisObj} is a function that calculates the return levels for a Generalized Pareto Distribution (GPD) using the parameters obtained from an analysis object. It takes into account non-stationarity by considering time-varying parameters and their associated standard errors.
#'
#' @param nonStationaryEvaParams The non-stationary parameters obtained from the analysis object.
#' @param returnPeriodsInYears The return periods for which to compute the return levels, expressed in years.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list with the following components:
#'   \itemize{
#'     \item \code{returnLevels} A vector of return levels corresponding to the specified return periods.
#'     \item \code{returnLevelsErrFit} A vector of standard errors for the return levels estimated based on the fit.
#'     \item \code{returnLevelsErrTransf} A vector of standard errors for the return levels estimated based on the transformed parameters.
#'   }
#'
#' @seealso \code{\link{tsEvaComputeReturnLevelsGPD}}, \code{\link{tsEasyParseNamedArgs}}
#'
#'
#' @family MyPackage functions
#'
#' @docType methods
#' @name tsEvaComputeReturnLevelsGPDFromAnalysisObj
#' @rdname tsEvaComputeReturnLevelsGPDFromAnalysisObj
#' @aliases tsEvaComputeReturnLevelsGPDFromAnalysisObj
tsEvaComputeReturnLevelsGPDFromAnalysisObj <- function(nonStationaryEvaParams, returnPeriodsInYears, ...) {
  args <- list(timeIndex = -1)
  varargin <- list(...)
  args <- tsEasyParseNamedArgs(varargin, args)
  timeIndex <- args$timeIndex
  epsilon <- nonStationaryEvaParams[[2]]$parameters$epsilon
  epsilonStdErr <- nonStationaryEvaParams[[2]]$paramErr$epsilonErr
  epsilonStdErrFit <- epsilonStdErr
  epsilonStdErrTransf <- 0
  thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
  thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
  timeHorizonInYears <- (thEnd - thStart) / 365.2425
  nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
  nonStationary <- "sigmaErrTransf" %in% names(nonStationaryEvaParams[2]$paramErr)
  if (timeIndex > 0) {
    sigma <- nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex]
    threshold <- nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex]
    sigmaStdErr <- nonStationaryEvaParams[[2]]$paramErr$sigmaErr[timeIndex]
    sigmaStdErrFit <- nonStationaryEvaParams[[2]]$paramErr$sigmaErrFit[timeIndex]
    sigmaStdErrTransf <- nonStationaryEvaParams[[2]]$paramErr$sigmaErrTransf[timeIndex]
    thresholdStdErr <- nonStationaryEvaParams[[2]]$paramErr$thresholdErr[timeIndex]
    thresholdStdErrFit <- 0
    thresholdStdErrTransf <- nonStationaryEvaParams[[2]]$paramErr$thresholdErrTransf[timeIndex]
  } else {
    sigma <- nonStationaryEvaParams[[2]]$parameters$sigma
    threshold <- nonStationaryEvaParams[[2]]$parameters$threshold
    sigmaStdErr <- nonStationaryEvaParams[[2]]$paramErr$sigmaErr
    if (nonStationary) {
      thresholdStdErr <- nonStationaryEvaParams[[2]]$paramErr$thresholdErr
      sigmaStdErrFit <- nonStationaryEvaParams[[2]]$paramErr$sigmaErrFit
      sigmaStdErrTransf <- nonStationaryEvaParams[[2]]$paramErr$sigmaErrTransf
      thresholdStdErrFit <- nonStationaryEvaParams[[2]]$paramErr$thresholdErrFit
      thresholdStdErrTransf <- nonStationaryEvaParams[[2]]$paramErr$thresholdErrTransf
    } else {
      thresholdStdErr <- 0
    }
  }
  returnLevels <- tsEvaComputeReturnLevelsGPD(epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, timeHorizonInYears, returnPeriodsInYears)
  if (nonStationary) {
    returnLevelsErrFit <- tsEvaComputeReturnLevelsGPD(epsilon, sigma, threshold, epsilonStdErrFit, sigmaStdErrFit, rep(0, length(thresholdStdErr)), nPeaks, timeHorizonInYears, returnPeriodsInYears)
    returnLevelsErrTransf <- tsEvaComputeReturnLevelsGPD(epsilon, sigma, threshold, epsilonStdErrTransf, sigmaStdErrTransf, thresholdStdErrTransf, nPeaks, timeHorizonInYears, returnPeriodsInYears)
  } else {
    returnLevelsErrFit <- returnLevelsErr
    returnLevelsErrTransf <- rep(0, length(returnLevelsErr))
  }
}


#' tsEvaSampleData Function
#'
#' This function calculates various statistics and data for time series evaluation.
#'
#' @param ms A matrix containing the time series data.
#' @param meanEventsPerYear The mean number of events per year.
#' @param minEventsPerYear The minimum number of events per year.
#' @param minPeakDistanceInDays The minimum peak distance in days.
#' @param tail The tail to be studied for POT selection, either 'high' or 'low'.
#'
#' @return A list containing the following elements:
#'   \itemize{
#'     \item \code{completeSeries}: The complete time series data.
#'     \item \code{POT}: The data for Peaks Over Threshold (POT) analysis.
#'     \item \code{years}: The years in the time series data.
#'     \item \code{Percentiles}: The desired percentiles and their corresponding values.
#'     \item \code{annualMax}: The annual maximum values.
#'     \item \code{annualMaxDate}: The dates corresponding to the annual maximum values.
#'     \item \code{annualMaxIndx}: The indices of the annual maximum values.
#'     \item \code{monthlyMax}: The monthly maximum values.
#'     \item \code{monthlyMaxDate}: The dates corresponding to the monthly maximum values.
#'     \item \code{monthlyMaxIndx}: The indices of the monthly maximum values.
#'   }
#' @seealso [tsGetPot()]
#' @import stats
#' @examples
#' # Generate sample data
#' data <- matrix(c(1:100, 101:200), ncol = 2)
#' colnames(data) <- c("Date", "Value")
#' # Calculate statistics and data
#' result <- tsEvaSampleData(data, 10, 5, 7, "high")
#'
#' # View the result
#' print(result)
#'
#' @export
tsEvaSampleData <- function(ms, meanEventsPerYear,minEventsPerYear, minPeakDistanceInDays,tail=NA) {

  pctsDesired = c(90, 95, 99, 99.9)

  args <- list(meanEventsPerYear = meanEventsPerYear,
               minEventsPerYear = minEventsPerYear,
               potPercentiles = c(seq(70,90,by=1), seq(91,99.5,by=0.5)))
  meanEventsPerYear = args$meanEventsPerYear
  minEventsPerYear = args$minEventsPerYear
  potPercentiles = args$potPercentiles
  if(is.na(tail)) print("haz for POT selection needs to be 'high' or 'low'")

  POTData <- tsGetPOT(ms, potPercentiles, meanEventsPerYear,minEventsPerYear,minPeakDistanceInDays, tail)

  vals <- quantile(ms[,2], pctsDesired/100,na.rm=T)
  percentiles <- list(precentiles = pctsDesired, values = vals)

  pointData <- list()
  pointData$completeSeries <- ms
  pointData$POT <- POTData

  pointDataA <- computeAnnualMaxima(ms)
  pointDataM <- computeMonthlyMaxima(ms)

  yrs <- unique(as.numeric(format(as.Date(ms[,1]+3600), "%Y")))
  yrs <- yrs - min(yrs)
  pointData$years <- seq(min(yrs),max(yrs),1)

  pointData$Percentiles <- percentiles
  pointData$annualMax=pointDataA$annualMax
  pointData$annualMaxDate=pointDataA$annualMaxDate
  pointData$annualMaxIndx=pointDataA$annualMaxIndx
  pointData$monthlyMax=pointDataM$monthlyMax
  pointData$monthlyMaxDate=pointDataM$monthlyMaxDate
  pointData$monthlyMaxIndx=pointDataM$monthlyMaxIndx

  return(pointData)
}



#' tsGetPOT Function
#'
#' This function calculates the Peaks Over Threshold (POT) for a given time series data.
#'
#' @param ms A matrix containing the time series data with two columns: the first column represents the time and the second column represents the values.
#' @param pcts A numeric vector specifying the percentiles to be used as thresholds for identifying peaks.
#' @param desiredEventsPerYear The desired number of events per year.
#' @param minEventsPerYear The minimum number of events per year.
#' @param minPeakDistanceInDays The minimum distance between two peaks in days.
#' @param tail The tail to be studied for POT selection, either 'high' or 'low'.
#'
#' @return A list containing the following fields:
#' \describe{
#'   \item{threshold}{The threshold value used for identifying peaks.}
#'   \item{thresholdError}{The error associated with the threshold value.}
#'   \item{percentile}{The percentile value used as the threshold.}
#'   \item{peaks}{The values of the identified peaks.}
#'   \item{stpeaks}{The start indices of the identified peaks.}
#'   \item{endpeaks}{The end indices of the identified peaks.}
#'   \item{ipeaks}{The indices of the identified peaks.}
#'   \item{time}{The time values corresponding to the identified peaks.}
#'   \item{pars}{The parameters of the Generalized Pareto Distribution (GPD) fitted to the peaks.}
#' }
#'
#' @import stats
#' @importFrom POT fitgpd
#' @seealso [tsSampleData()]
#' @examples
#' # Create a sample time series data
#' time <- seq(as.Date("2000-01-01"), as.Date("2000-12-31"), by = "day")
#' values <- rnorm(length(time))
#' ms <- cbind(time, values)
#'
#' # Calculate the POT using the tsGetPOT function
#' pcts <- c(90, 95, 99)
#' desiredEventsPerYear <- 5
#' minEventsPerYear <- 2
#' minPeakDistanceInDays <- 10
#' mode <- "flood"
#' POTdata <- tsGetPOT(ms, pcts, desiredEventsPerYear, minEventsPerYear, minPeakDistanceInDays, mode)
#'
#' # Print the results
#' print(POTdata)
#'
#' @export
tsGetPOT <- function(ms, pcts, desiredEventsPerYear,minEventsPerYear, minPeakDistanceInDays, tail) {

  if (minPeakDistanceInDays == -1) {
    stop("label parameter 'minPeakDistanceInDays' must be set")
  }
  dt1=min(diff(ms[,1]),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  if (tdim=="hours") dt=dt/24
  if (tdim=="seconds") dt=dt/3600
  minPeakDistance <- minPeakDistanceInDays/dt
  minRunDistance <- minPeakDistance
  nyears <- round(as.numeric((max(ms[,1]) - min(ms[,1]))/365.25))
  if (length(pcts) == 1) {
    pcts = c(pcts - 3, pcts)
    desiredEventsPerYear = -1
  }

  numperyear <- rep(NA, length(pcts))
  minnumperyear <- rep(NA, length(pcts))
  thrsdts <- rep(NA, length(pcts))
  gpp=rep(NA, length(pcts))
  devpp=rep(NA, length(pcts))
  dej=0
  skip=0
  trip=NA
  perfpen=0
  for (ipp in 1:length(pcts)) {
    #Skip is used to prevent finding peaks for unappropriate thresholds
    if (skip>0) {
      skip=skip-1
    }else{
      if(dej==0){
        thrsdt <- quantile(ms[,2],pcts[ipp]/100,na.rm=T)
        thrsdts[ipp] <- thrsdt
        ms[,2][which(is.na(ms[,2]))]=-9999
        minEventsPerYear=1

        if(tail=="high") {
          #boundaries of shape parameter
          shape_bnd=c(-0.5,1)
          pks <- pracma::findpeaks(ms[,2],minpeakdistance = minPeakDistance, minpeakheight = thrsdt)
        }
        if(tail=="low") {
          pks <- declustpeaks(data = ms[,2] ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=thrsdt)
          shape_bnd=c(-1.5,0)
        }
        numperyear[ipp] <- length(pks[,1])/nyears
        if(numperyear[ipp]>=3*desiredEventsPerYear & ipp<(length(pcts)-5)) skip = floor(length(pcts)/8)
        if(numperyear[ipp]<minEventsPerYear) {
          perfpen=(pcts[ipp]*100)
        }
        if(numperyear[ipp]<(0.7*minEventsPerYear)) {
          perfpen=9999*(pcts[ipp]*100)
        }
        if(numperyear[ipp]<=desiredEventsPerYear+1 & dej==0){
          fgpd=suppressWarnings(POT::fitgpd(pks[,1], threshold = thrsdt, est = "mle",method="BFGS",std.err.type = "expected"))
          gpdpar=fgpd$fitted.values
          deviance=fgpd$deviance
          devpp[ipp]=AIC(fgpd)+perfpen
          gpp[ipp]=gpdpar[2]
          nperYear <- tsGetNumberPerYear(ms, pks[,2])
          minnumperyear[ipp] <- min(nperYear$Freq, na.rm = TRUE)
        }
      }
    }
  }
  devpp[1]=NA
  print(shape_bnd)
  if(is.na(trip)){
    isok=F
    count=0
    devpx=devpp
    #discard the fits that provide negative values
    while(isok==F){
      trip=which.min(devpx)
      isok=between(gpp[trip], shape_bnd[1], shape_bnd[2])
      count=count+1
      if(isok==F)devpx[trip]=devpx[trip]+9999
      if(count>(length(devpx)-1)){
        #safety measure for stability of parameter
        dshap=c(0,diff(gpp))
        #Penalizing fits with positive shape parameters for low tail
        if(tail=="low") devpp[which(gpp>=0)]=devpp[which(gpp>=0)]+99999
        devpp[which(abs(dshap)>0.5)]=devpp[which(abs(dshap)>0.5)]+9999
        trip=which.min(devpp)
        print(paste0("shape outside boudaries: ",round(gpp[trip],2)))
        isok=T
      }
    }
  }
  cat(paste0("\nmax threshold is: ", pcts[trip],"%"))
  cat(paste0("\naverage number of events per year = ",round(numperyear[trip],1) ))

  diffNPerYear <- mean(diff(na.omit(rev(numperyear)), na.rm = TRUE))
  if (diffNPerYear == 0) diffNPerYear <- 1
  diffNPerYear <- 1
  thresholdError <- -mean(diff(na.omit(thrsdts))/diffNPerYear)/2
  indexp <- trip
  if (!is.na(indexp)) {
    thrsd <- quantile(ms[,2],pcts[indexp]/100)
    pct <- pcts[indexp]
  } else {
    thrsd <- 0
    pct
  }
  # Find peaks in the second column of the matrix 'ms'
  if(tail=="high") pks_and_locs <- pracma::findpeaks(ms[,2],minpeakdistance = minPeakDistance, minpeakheight = thrsd)
  if(tail=="low") pks_and_locs <- declustpeaks(data = ms[,2] ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=thrsd)

  # Assign peaks and peak locations to separate variables
  pks <- pks_and_locs[,1]
  locs <- pks_and_locs[,2]
  st<-pks_and_locs[,3]
  end=pks_and_locs[,4]

  # Create a list to store results
  POTdata <- list()
  # Assign values to the fields of the list
  POTdata[['threshold']] <- thrsd
  POTdata[['thresholdError']] <- thresholdError
  POTdata[['percentile']] <- pct
  POTdata[['peaks']] <- pks
  POTdata[['stpeaks']] <- st
  POTdata[['endpeaks']] <- end
  POTdata[['ipeaks']] <- locs
  POTdata[['time']] <- ms[locs, 1]
  POTdata[['pars']] <- gpdpar


  return(POTdata)
}

#' Calculate GEV and GPD statistics and return levels
#'
#' This function calculates the Generalized Extreme Value (GEV) and
#' Generalized Pareto Distribution (GPD) statistics and return levels
#' for a given dataset of extreme values.
#'
#' @param pointData A list containing the dataset of extreme values. It should include the following components:
#'   \describe{
#'     \item{annualMax}{A vector of annual maximum values.}
#'     \item{annualMaxDate}{A vector of dates corresponding to the annual maximum values.}
#'     \item{monthlyMax}{A matrix of monthly maximum values.}
#'   }
#' @param alphaCI The confidence level for the confidence intervals of the parameter estimates. Default is 0.95.
#' @param gevMaxima The type of maxima to use for GEV fitting. Can be either 'annual' or 'monthly'. Default is 'annual'.
#' @param gevType The type of GEV distribution to use. Can be either 'GEV', 'Gumbel'. Default is 'GEV'.
#' @param evdType The types of extreme value distributions to calculate. Can be a combination of 'GEV' and 'GPD'. Default is c('GEV', 'GPD').
#' @param shape_bnd The lower and upper bounds for the shape parameter of the GEV distribution. Default is c(-0.5, 1).
#'
#' @importFrom evd fgev
#' @return A list containing the following components:
#'   \describe{
#'     \item{EVmeta}{A list containing metadata about the analysis. It includes the following components:
#'       \describe{
#'         \item{Tr}{A vector of return periods for which return levels are calculated.}
#'       }
#'     }
#'     \item{EVdata}{A list containing the calculated statistics and return levels. It includes the following components:
#'       \describe{
#'         \item{GEVstat}{A list containing the GEV statistics and return levels. It includes the following components:
#'           \describe{
#'             \item{method}{The method used for fitting the GEV distribution.}
#'             \item{values}{A vector of return levels calculated using the GEV distribution.}
#'             \item{parameters}{A vector of parameter estimates for the GEV distribution.}
#'             \item{paramCIs}{A matrix of confidence intervals for the parameter estimates.}
#'           }
#'         }
#'         \item{GPDstat}{A list containing the GPD statistics and return levels. It includes the following components:
#'           \describe{
#'             \item{method}{The method used for fitting the GPD distribution.}
#'             \item{values}{A vector of return levels calculated using the GPD distribution.}
#'             \item{parameters}{A vector of parameter estimates for the GPD distribution.}
#'             \item{paramCIs}{A matrix of confidence intervals for the parameter estimates.}
#'           }
#'         }
#'       }
#'     }
#'     \item{isValid}{A logical value indicating whether the analysis was performed or not.}
#'   }
#'
#' @examples
#' # Create a sample dataset
#' pointData <- list(
#'   annualMax = c(10, 15, 12, 18, 20),
#'   annualMaxDate = c("2010-01-01", "2011-01-01", "2012-01-01", "2013-01-01", "2014-01-01"),
#'   monthlyMax = matrix(c(5, 8, 6, 10, 7, 12, 9, 14, 11, 16, 13, 18), ncol = 3)
#' )
#'
#' # Calculate GEV and GPD statistics
#' result <- tsEVstatistics(pointData)
#' result$EVdata$GEVstat$values
#' result$EVdata$GPDstat$values
#'
#' @export
tsEVstatistics <- function(pointData, alphaCI = 0.95, gevMaxima = 'annual', gevType = 'GEV', evdType = c('GEV', 'GPD'),shape_bnd=c(-0.5,1)) {
  # Create empty data structures
  EVmeta <- list()
  EVdata <- list()
  isValid <- TRUE

  minGEVSample <- 7

  if(is.null(alphaCI)){alphaCI <- .95}
  if(is.null(gevMaxima)){gevMaxima <- 'annual'}
  if(is.null(gevType)){gevType <- 'GEV'}
  if(is.null(evdType)){evdType <- c('GEV', 'GPD')}

  # Define Tr vector
  Tr <- c(5,10,20,50,100,200,500,1000)
  EVmeta$Tr <- Tr
  nyears <- length(pointData$annualMax)[1]

  imethod <- 1
  methodname <- 'GEVstat'
  paramEsts <- numeric(3)
  paramCIs <- matrix(NA, nrow = 2, ncol = 3)
  rlvls <- numeric(length(Tr))

  # GEV statistics
  if (('GEV' %in% evdType) && !is.null(pointData$annualMax)) {
    if (gevMaxima == 'annual') {
      tmpmat <- pointData$annualMax
    } else if (gevMaxima == 'monthly') {
      tmpmat <- pointData$monthlyMax
    } else {
      stop(paste0('Invalid gevMaxima type: ', gevMaxima))
    }
    iIN <- length(tmpmat)

    if (sum(iIN) >= minGEVSample) {
      tmp <- data.frame(yr=year(pointData$annualMaxDate),dt=tmpmat)
      # Perform GEV/Gumbel fitting and computation of return levels
      if (gevType == "GEV"){
        stdfit=TRUE
        #try to fit GEV with bounded shape parameters and stderr, reduces the constrains if no fit.
        fit <- suppressWarnings(try(evd::fgev(x=tmp$dt,method="L-BFGS-B",lower=c(-Inf,-Inf,shape_bnd[1]),upper=c(Inf,Inf,shape_bnd[2]),std.err = T),TRUE))
        if(inherits(fit, "try-error")){
          stdfit=FALSE
          print("Not able to fit GEV stderr")
          fit <- suppressWarnings(try(evd::fgev(x=tmp$dt,method="L-BFGS-B",lower=c(-Inf,-Inf,shape_bnd[1]),upper=c(Inf,Inf,shape_bnd[2]),std.err = F),TRUE))
        }
        if(inherits(fit, "try-error")){
          stdfit=FALSE
          print("Not able to fit GEV with constrained parameters")
          fit <- suppressWarnings(try(evd::fgev(x=tmp$dt,std.err = F),TRUE))
        }
        paramEsts <- c(mu=fit$par[1],sigma=fit$par[2],xi=fit$par[3])
        alphaCIx=1-alphaCI
        if (stdfit==TRUE){
          probs <- c(alphaCIx/2, 1-alphaCIx/2)
          # Compute the CI for k using a normal distribution for khat.
          kci <- try(qnorm(probs, paramEsts[3], fit$std.err[3]),TRUE)
          kci[kci < -1] <- -1
          # Compute the CI for sigma using a normal approximation for log(sigmahat)
          # and transform back to the original scale.
          lnsigci <- try(qnorm(probs, log(paramEsts[2]), fit$std.err[2]/paramEsts[2]),silent=T)
          muci <- qnorm(probs, paramEsts[1], fit$std.err[1])
          paramCIs <- cbind(muci=muci, sigci=exp(lnsigci), kci=kci)
        }else{
          #not able to generate parameters CI
          paramCIs <- cbind(muci=NA, sigci=NA, kci=NA)
        }
      } else if (gevType == "gumbel"){
        fit <- texmex::evm(y = dt, data = tmp, family = gumbel)
        paramEsts <- c(fit$par[1], exp(fit$par[2]),0)

        alphaCIx=1-alphaCI
        probs <- c(alphaCIx/2, 1-alphaCIx/2)
        # Compute the CI for k using a normal distribution for khat.
        if(is.character(fit$se[1])){
          print("Gumbel fitted")
          #methodname <- 'No fit'
          gevType = "gumbel"
        }else{
          print("Gumbel fitted")
          #gevType = "gumbel"
          kci <- c(NA,NA)
          # Compute the CI for sigma using a normal approximation for log(sigmahat)
          # and transform back to the original scale.
          lnsigci <- try(qnorm(probs, log(paramEsts[2]), fit$se[2]))

          muci <- qnorm(probs, paramEsts[1], fit$se[1])
          paramCIs <- cbind(muci=muci, sigci=exp(lnsigci), kci=kci)
        }
      }
      else {
        stop("tsEVstatistics: invalid gevType: ", gevType, ". Can be only GEV, or Gumbel")
      }
      rlvls <- qgev(1-1/Tr, paramEsts[1], paramEsts[2], paramEsts[3])

    }else{
      print("Could not fit GEV")
      methodname <- 'No fit'
      rlvls=NA
      paramEsts=NA
      paramCIs=NA

    }
  }

  EVdata$GEVstat <- list(method=methodname,
                         values=rlvls,
                         parameters=paramEsts,
                         paramCIs = paramCIs)

  # Create output structures for GEV statistics


  # GPD statistics
  imethod <- 2
  methodname <- 'GPDstat'
  paramEsts <- numeric(6)
  paramCIs <- matrix(NA, nrow = 2, ncol = 3)
  rlvls <- numeric(length(Tr))

  if (('GPD' %in% evdType) && !is.null(pointData$annualMax)) {
    # Perform GPD fitting and computation of return levels
    print("Fitted GPD")
    ik <- 1
    th=pointData$POT$threshold
    d1 <- pointData$POT$peaks


    fit <- suppressWarnings(try(POT::fitgpd(d1, threshold = th, est = "mle",
                                       method="BFGS",std.err.type = "expected"),TRUE))
    if(!inherits(fit, "try-error")){
      ksi <- fit$par[2]
      sgm <- fit$par[1]
      alphaCIx=1-alphaCI
      probs <- c(alphaCIx/2, 1-alphaCIx/2)
      kci <- try(qnorm(probs, ksi, fit$std.err[2]),silent=T)
      kci[kci < -1] <- -1
      # Compute the CI for sigma using a normal approximation for log(sigmahat)
      # and transform back to the original scale.
      lnsigci <- try(qnorm(probs, log(sgm), fit$std.err[1]/sgm))
      paramCIs <- cbind(kci, sigci=exp(lnsigci))

      # Create output structures for GPD statistics

      # Assign values to paramEstsall
      paramEstsall <- c(sgm, ksi, pointData$POT$threshold, length(d1), length(pointData$POT$peaks), pointData$POT$percentile)
      # Assign values to rlvls
      rlvls <- pointData$POT$threshold + (sgm/ksi) * ((((length(d1)/length(pointData$POT$peaks))*(1/Tr))^(-ksi))-1)

      EVdata$GPDstat <- list(method=methodname,
                             values=rlvls,
                             parameters=paramEstsall,
                             paramCIs = paramCIs)
    }else {
      methodname <- 'No fit'
      ik <- 1
      th=pointData$POT$threshold
      d1 <- pointData$POT$peaks
      paramEstsall <- c(pointData$POT$pars[1], pointData$POT$pars[2],
                        pointData$POT$threshold, length(d1),
                        length(pointData$POT$peaks), pointData$POT$percentile)
      EVdata$GPDstat <-  list(method=methodname,
                              values=NA,
                              parameters=paramEstsall,
                              paramCIs = NA)
      print("could not estimate GPD: bounded distribution")
    }
  } else {
    methodname <- 'No fit'
    ik <- 1
    th=pointData$POT$threshold
    d1 <- pointData$POT$peaks
    paramEstsall <- c(pointData$POT$pars[1], pointData$POT$pars[2],
                      pointData$POT$threshold, length(d1),
                      length(pointData$POT$peaks), pointData$POT$percentile)
    EVdata$GPDstat <-  list(method=methodname,
                            values=NA,
                            parameters=paramEstsall,
                            paramCIs = NA)
    print("could not estimate GPD: bounded distribution")
  }

  # Return outputs
  return(list(EVmeta = EVmeta, EVdata = EVdata, isValid = isValid))
}



# Helpers for main functions----------------

#' Get the number of events per year
#'
#' This function calculates the number of events per year based on a given time series and a set of locations.
#'
#' @param ms A data frame representing the time series data, where the first column contains the dates of the events.
#' @param locs A vector of indices representing the locations of interest in the time series.
#'
#' @return A data frame with two columns: "year" and "Freq". The "year" column contains the years, and the "Freq" column contains the number of events per year.
#'
#' @examples
#' # Create a sample time series data frame
#' ms <- data.frame(date = seq(as.Date("2000-01-01"), as.Date("2022-12-31"), by = "day"))
#'
#' # Generate random events
#' set.seed(123)
#' events <- sample(ms$date, 1000)
#'
#' # Get the number of events per year
#' tsGetNumberPerYear(ms, events)
#'
#' @importFrom dplyr full_join
#' @importFrom lubridate year
#' @export
tsGetNumberPerYear <- function(ms, locs){

  # Get all years
  sdfull <- seq(min(ms[,1]), max(ms[,1]), by = "days")

  # Make full time vector
  sdfull2 <- sdfull - min(sdfull)

  # Round to years
  sdfull2 <- year(sdfull)
  sdfull2p=data.frame(year=unique(sdfull2))


  # Prepare the series time vector
  sdser <- ms[,1] - min(ms[,1])
  sdser <- year(ms[,1])

  # Get years vector
  sdp <- ms[locs,1]
  sdp <- year(as.Date(sdp))


  # Get number of events per year
  nperYear <- data.frame(table(sdp))
  nperYear$year=as.numeric(as.character(nperYear$sdp))
  nperYear=dplyr::full_join(nperYear,sdfull2p,by=c("year"="year"))
  shit=table(sdser)
  nperYear$Freq[which(is.na(nperYear$sdp))]=0
  nperYear$opy <- table(sdser)
  nperYear$Freq[which(nperYear$opy<184)] <- 0
  return(nperYear)
}


#' Find the index of the yearly maximum value in a subset of a vector
#'
#' This function takes a subset of a vector and returns the index of the maximum value in that subset.
#'
#' @param subIndxs A numeric vector representing the subset of indices to consider.
#' @return The index of the maximum value in the subset.
#' @examples
#' srs <- c(10, 20, 30, 40, 50)
#' findYMax(c(1, 3, 5))
#' # Output: 5
findYMax <- function(subIndxs) {
  subIndxMaxIndx <- which.max(srs[subIndxs])
  subIndxs[subIndxMaxIndx]
}

#' Find the index of the monthly maximum value in a given subset of a vector.
#'
#' This function takes a vector and a subset of indices and returns the index of the maximum value within that subset.
#'
#' @param indxs A numeric vector specifying the subset of indices.
#' @return The index of the maximum value within the specified subset.
#' @examples
#' srs <- c(10, 20, 30, 40, 50)
#' findMMax(c(1, 3, 5))
#' # Output: 5
findMMax <- function(indxs) {
  maxVal <- max(srs[indxs])
  maxIndx <- indxs[which.max(srs[indxs])]
  return(maxIndx)
}

#' Compute Annual Maxima
#'
#' This function computes the annual maxima of a time series.
#'
#' @param timeAndSeries A matrix or data frame containing the time stamps and series values.
#'                      The first column should contain the time stamps and the second column should contain the series values.
#'
#' @return A list containing the annual maximum values, their corresponding dates, and their indices.
#'         - \code{annualMax}: A numeric vector of annual maximum values.
#'         - \code{annualMaxDate}: A vector of dates corresponding to the annual maximum values.
#'         - \code{annualMaxIndx}: A vector of indices indicating the positions of the annual maximum values in the original time series.
#'
#' @examples
#' timeAndSeries <- data.frame(timeStamps = c("2021-01-01", "2021-01-02", "2021-01-03"),
#'                             series = c(10, 15, 8))
#' computeAnnualMaxima(timeAndSeries)
#'
#' @export
computeAnnualMaxima <- function(timeAndSeries) {
  timeStamps <- timeAndSeries[,1]
  dt1=min(diff(timeStamps),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  tmvec <- as.Date(timeStamps)
  if (tdim!="days")   tmvec <- as.Date(timeStamps+3600)
  srs <- timeAndSeries[,2]
  years <- year(tmvec)

  annualMaxIndx <- tapply(1:length(srs), years, findMax)
  annualMaxInx <- annualMaxIndx[!is.na(annualMaxIndx)]
  annualMaxIndx=as.vector(unlist((annualMaxIndx)))
  annualMax <- srs[annualMaxIndx]
  annualMaxDate <- timeStamps[annualMaxIndx]

  return(list(annualMax = annualMax, annualMaxDate = annualMaxDate, annualMaxIndx = annualMaxIndx))
}




#' Compute Monthly Maxima
#'
#' This function computes the monthly maxima of a time series.
#'
#' @param timeAndSeries A data frame containing the time stamps and series values.
#'                      The first column should contain the time stamps, and the second column should contain the series values.
#'
#' @return A list containing the monthly maxima, corresponding dates, and indices.
#'         - \code{monthlyMax}: A vector of the monthly maximum values.
#'         - \code{monthlyMaxDate}: A vector of the dates corresponding to the monthly maximum values.
#'         - \code{monthlyMaxIndx}: A vector of the indices of the monthly maximum values in the original series.
#'
#' @examples
#' timeAndSeries <- data.frame(timeStamps = c("2022-01-01", "2022-01-02", "2022-02-01", "2022-02-02"),
#'                             series = c(10, 15, 5, 20))
#' computeMonthlyMaxima(timeAndSeries)
#'
#' @export
computeMonthlyMaxima<- function(timeAndSeries) {
  timeStamps <- timeAndSeries[,1]
  srs <- timeAndSeries[,2]

  tmvec <- as.Date(timeStamps+3600)
  yrs <- as.integer(format(tmvec, "%Y"))
  mnts <- as.integer(format(tmvec, "%m"))
  mnttmvec <- data.frame(yrs, mnts)
  valsIndxs <- 1:length(srs)

  monthlyMaxIndx <- aggregate(valsIndxs ~ yrs + mnts, data = mnttmvec, findMax)
  monthlyMaxIndx$valsIndxs=as.numeric(monthlyMaxIndx$valsIndxs)
  monthlyMaxIndx <- monthlyMaxIndx[order(monthlyMaxIndx$yrs, monthlyMaxIndx$mnts), "valsIndxs"]
  monthlyMax <- srs[monthlyMaxIndx]
  monthlyMaxDate <- timeStamps[monthlyMaxIndx]

  return(list(monthlyMax = monthlyMax, monthlyMaxDate = monthlyMaxDate, monthlyMaxIndx = monthlyMaxIndx))
}


#' declustpeaks Function
#'
#' This function takes in a data vector, minimum peak distance, minimum run distance, and a threshold value.
#' It finds peaks in the data vector based on the minimum peak distance and threshold value.
#' It then declusters the peaks based on the minimum run distance and threshold value.
#' The function returns a data frame with information about the peaks, including the peak value,
#' start and end indices, duration, and cluster information.
#'
#' @param data A numeric vector representing the data.
#' @param minpeakdistance An integer specifying the minimum distance between peaks.
#' @param minrundistance An integer specifying the minimum distance between runs.
#' @param qt A numeric value representing the threshold for peak detection.
#'
#' @import texmex stats
#' @return A data frame with information about the peaks, including the peak value,
#' start and end indices, duration, and cluster information.
#'
#' @examples
#' data <- c(1, 2, 3, 4, 5, 4, 3, 2, 1)
#' declustpeaks(data, minpeakdistance = 2, minrundistance = 2, qt = 3)
#'
#' @export
declustpeaks<-function(data,minpeakdistance = 10 ,minrundistance = 7, qt){
  pks <- pracma::findpeaks(data,minpeakdistance = minpeakdistance, minpeakheight = qt)
  peakev=texmex::declust(data,threshold=qt, r=minrundistance)

  Qval=peakev$thExceedances
  intcl=c(TRUE,peakev$InterCluster)
  peakex=data.frame(Qval,intcl, peakev$clusters,peakev$isClusterMax,peakev$exceedanceTimes)
  names(peakex)=c("Qs","Istart","clusters","IsClustermax","exceedances")
  evmax=peakex[which(peakex$IsClustermax==T),]
  ziz=aggregate(peakex$exceedance,
                by=list(clust=peakex$clusters), FUN=function(x) c(max(x)-min(x)+1))
  st=stats::aggregate(peakex$exceedance,
               by=list(clust=peakex$clusters), FUN=function(x) c(min(x)))
  end=stats::aggregate(peakex$exceedance,
                by=list(clust=peakex$clusters), FUN=function(x) c(max(x)))
  evmax$dur=ziz$x
  evmax$durx=peakev$sizes
  evmax$stdate=peakex$date[which(evmax$Istart==T)]
  evmax$cm = peakex$exceedanceTimes[peakev$isClusterMax]
  evmax$Qv = peakex$thExceedances[peakev$isClusterMax]
  peakdt=data.frame(evmax$Qs,evmax$exceedances,st$x,end$x,ziz$x,evmax$clusters)
  names(peakdt)=c("Q","max","start","end","dur","cluster")
  peakdt=peakdt[order(peakdt$Q,decreasing = TRUE),]

  return(peakdt)
}






# empirical probability for return level plots
#' Empirical Distribution Function
#'
#' This function calculates the empirical distribution function for a given dataset.
#'
#' @param x A numeric vector representing the dataset.
#' @param nyr An integer representing the number of years in the dataset.
#'
#' @return A data frame containing the empirical return periods, hazard return periods,
#' Cunnane return periods, Gumbel values, empirical probabilities, hazard probabilities,
#' Cunnane probabilities, the original dataset, and the timestamp.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' nyr <- 5
#' empdis(x, nyr)
#'
#' @export
empdis <- function(x, nyr) {
  ts <- seq(1:length(x))
  ts <- as.data.frame(ts[order(x)])
  x <- as.data.frame(x[order(x)])
  epyp <- length(x[, 1]) / nyr
  rank <- apply(x, 2, rank, na.last = "keep")
  hazen <- (rank - 0.5) / length(x[, 1])
  cunnane <- (rank - 0.4) / (length(x[, 1]) + 0.2)
  gumbel <- -log(-log(hazen))
  n <- length(x$`x[order(x)]`)
  rpc <- (1 / (1 - (1:n) / (n + 1))) / 12
  nasmp <- apply(x, 2, function(x) sum(!is.na(x)))
  epdp <- rank / rep(nasmp + 1, each = nrow(rank))
  empip <- data.frame(1 / (epyp * (1 - epdp)), 1 / (epyp * (1 - hazen)), 1 / (epyp * (1 - cunnane)), gumbel, epdp, hazen, cunnane, x, ts)
  names(empip) <- c("emp.RP", "haz.RP", "cun.RP", "gumbel", "emp.f", "emp.hazen", "emp.cunnane", "Q", "timestamp")
  return(empip)
}


#' Empirical Distribution Function
#'
#' This function calculates the empirical distribution function for a given dataset,
#'  with a focus on low values
#'
#' @param x A numeric vector or matrix representing the discharge values.
#' @param nyr An integer representing the number of years.
#'
#' @return A data frame containing the following columns:
#'   \describe{
#'     \item{emp.RP}{Empirical return period.}
#'     \item{haz.RP}{Hazen return period.}
#'     \item{gumbel}{Gumbel frequency.}
#'     \item{emp.f}{Empirical frequency.}
#'     \item{emp.hazen}{Empirical Hazen frequency.}
#'     \item{Q}{Discharge values.}
#'   }
#'
#' @examples
#' x <- c(10, 20, 30, 40, 50)
#' nyr <- 5
#' empdisl(x, nyr)
#'
#' @export
empdisl <- function(x, nyr) {
  x <- as.data.frame(x[order(x, decreasing = TRUE)])
  epyp <- length(x[, 1]) / nyr
  rank <- apply(-x, 2, rank, na.last = "keep")
  hazen <- (rank - 0.5) / length(x[, 1])
  gumbel <- -log(-log(hazen))
  nasmp <- apply(x, 2, function(x) sum(!is.na(x)))
  epdp <- rank / rep(nasmp + 1, each = nrow(rank))
  empip <- data.frame(1 / (epyp * (1 - epdp)), 1 / (epyp * (1 - hazen)), gumbel, epdp, hazen, x)
  names(empip) <- c("emp.RP", "haz.RP", "gumbel", "emp.f", "emp.hazen", "Q")
  return(empip)
}
