# suppressWarnings(suppressMessages(library(tsibble)))
# suppressWarnings(suppressMessages(library(ggplot2)))
# suppressWarnings(suppressMessages(library(scales)))
# suppressWarnings(suppressMessages(library(pracma)))
# suppressWarnings(suppressMessages(library(lubridate)))
# suppressWarnings(suppressMessages(library(texmex)))
# suppressWarnings(suppressMessages(library(ggplot2)))
# suppressWarnings(suppressMessages(library(moments)))
# suppressWarnings(suppressMessages(library(dplyr)))
# suppressWarnings(suppressMessages(library(changepoint)))
# suppressWarnings(suppressMessages(library(iemisc)))
# suppressWarnings(suppressMessages(library(evd)))
# suppressWarnings(suppressMessages(library(POT)))

#Main function for non-stationary analysis
#' TsEvaNs Function
#'
#' This function performs non-stationary extreme value analysis (EVA) on a time series data.
#'
#' @param timeAndSeries A data frame containing the timestamp and corresponding series data.
#' @param timeWindow The time window for analysis.
#' @param ... Additional arguments for customization.
#' @param transfType The transformation type for non-stationary EVA. It can be one of the following:
#' \itemize{trend: Long-term variations of the timeseries.
#' seasonal: Long-term and seasonal variations of extremes.
#' trendCIPercentile: Long-term variations of extremes using a specified percentile.
#' trendPeak: Long-term variations of the peaks.}
#' @param minPeakDistanceInDays The minimum peak distance in days.
#' @param seasonalityVar A logical value indicating whether to consider seasonality in the analysis.
#' @param potEventsPerYear The number of potential events per year.
#' @param minEventsPerYear The minimum number of events per year.
#' @param gevMaxima The type of maxima for the GEV distribution (annual or monthly).
#' @param ciPercentile The percentile value for confidence intervals.
#' @param gevType The type of GEV distribution (GEV or GPD).
#' @param evdType The type of extreme value distribution (GEV or GPD).
#' @param tail The mode of the analysis (e.g., high for flood peaks or low for drought peaks).
#' @param epy The average number of events per year, can be specified by the
#' user or automatically set according to the tail selected.
#' @param lowdt The temporal resolultion used for low values. default is 7 days.
#'
#' @return A list containing the results of the non-stationary EVA. Containing the following components:
#' \itemize{nonStationaryEvaParams: The estimated parameters for non-stationary EVA.
#' Parameters include GEV and GPD parameters for each timestep, confidence intervals, and other statistical measures.
#' stationaryTransformData: The transformed data for stationary EVA.Includes the stationary series, trend, and standard deviation series.}
#'
#' @details The function takes a time series data and performs non-stationary EVA using various transformation types and parameters depending on the input data provided.
#' Results are returned as a list containing the non-stationary EVA parameters and the transformed data for stationary EVA and can be used as input for further analysis.
#' In particular for the following function
#'
#' @examples
#' # Example usage of TsEvaNs function
#' TimeAndSeries <- data.frame(timestamp = seq(as.POSIXct("2000-01-01"), as.POSIXct("2020-12-31"), by = "day",data = rnorm(7671))
#' timeWindow <- 365
#' result <- TsEvaNs(timeAndSeries, timeWindow, transfType = 'trendPeaks', ciPercentile = 90, mode = 'drought')
#'
#' @export

TsEvaNs<- function(timeAndSeries, timeWindow, ...){
  ota=list(...)
  args <- list(transfType = 'trendCIPercentile',
               minPeakDistanceInDays = -1,
               seasonalityVar = F,
               potEventsPerYear = -1,
               minEventsPerYear = -1,
               gevMaxima = 'annual',
               ciPercentile = 90,
               gevType = 'GEV',
               evdType = c('GEV', 'GPD'),
               tail="low",
               epy=-1,
               trans=NULL,
               lowdt=7,
               TrendTh=0.75)

  args <- tsEasyParseNamedArgs(ota, args)
  seasonalityVar<-args$seasonalityVar
  ciPercentile <- args$ciPercentile
  transfType <- args$transfType
  evdType <- args$evdType
  gevType <- args$gevType
  TrendTh <- args$TrendTh
  epy <- args$epy
  tail=args$tail
  lowdt=args$lowdt
  minPeakDistanceInDays<-args$minPeakDistanceInDays
  trans=args$trans
  mode=args$mode


  timeStamps=as.POSIXct(timeAndSeries$timestamp)
  dt1=min(diff(timeStamps),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  if (tdim=="hours") dt=dt/24
  if (tdim=="seconds") dt=dt/3600
  series=timeAndSeries$data

  if (epy==-1){
    if (tail=="high") epy=3
    if (tail=="low") epy=2
  }

  #If the biggest event is more than 100 times greater than the second biggest event
  sanitycheck = computeAnnualMaxima(timeAndSeries);
  anmax=sanitycheck$annualMax[order(sanitycheck$annualMax,decreasing = T)]
  aloc=sanitycheck$annualMaxIndx[order(sanitycheck$annualMax,decreasing = T)][1]
  yrs=year(sanitycheck$annualMaxDate[order(sanitycheck$annualMax,decreasing = T)][1])
  x = anmax[1]/anmax[2]
  if (x > 100)
  {
    print(paste0("biggest event ", x, " times bigger than second biggest"))
    print ("removing this event from timeserie, reruning first steps")
    series[(aloc[1]-50):(aloc[1]+50)]=mean(series)
  }

  if ( transfType != 'trend' | transfType != 'seasonal' | transfType != 'trendCIPercentile' | transfType != 'seasonalCIPercentile' | transfType != 'trendPeak'){
    cat('\nnonStationaryEvaJRCApproach: transfType can be in (trend, seasonal, trendCIPercentile, trendPeak)')}

  if (minPeakDistanceInDays == -1) print('label parameter minPeakDistanceInDays must be set')

  nonStationaryEvaParams = c()
  stationaryTransformData = c()

  # default shape parameter bounds
  shape_bnd=c(-0.5,1)

  if (tail=="low"){
    #default 7-day flow for low flow, can be modified by user
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(series), by = lowdt/dt)
    series=series[indices_to_extract]
    timeStamps=timeStamps[indices_to_extract]
    shape_bnd=c(-1,0)
    if (trans=="rev"){
      series=-1*series
    }else if(trans=="inv"){
      series=1/series
    }else if (trans=="lninv"){
      series=-log(series)
    }

  }
  if (transfType == 'trend'){
    cat('\nevaluating long term variations of extremes')
    trasfData = tsEvaTransformSeriesToStationaryTrendOnly(timeStamps, series, timeWindow)
    gevMaxima = 'annual';
    potEventsPerYear = epy;
    minEventsPerYear = 1


  }else  if (transfType == 'trendChange'){
    cat('\nevaluating long term variations of extremes and change point detection')
    trasfData = tsEvaTransformSeriesToStationaryTrendAndChangepts(timeStamps, series, timeWindow)
    gevMaxima = 'annual';
    potEventsPerYear = epy;
    minEventsPerYear = 1

  }else if (transfType == 'seasonal'){
    cat('\nevaluating long term an seasonal variations of extremes')
    trasfData = tsEvaTransformSeriesToStationaryMultiplicativeSeasonality(timeStamps, series, timeWindow, seasonalityVar=seasonalityVar)
    gevMaxima = 'monthly'
    potEventsPerYear = 12
    minEventsPerYear = 12

  } else if (transfType == 'trendCIPercentile') {
    if (is.na(ciPercentile)){
      print('For trendCIPercentile transformation the label parameter cipercentile is mandatory')
    }
    cat(paste0('\nevaluating long term variations of extremes using the ', ciPercentile, 'th percentile'))
    trasfData = tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile( timeStamps, series, timeWindow, ciPercentile);
    gevMaxima = 'annual'
    potEventsPerYear = epy
    minEventsPerYear = 1

  }else if (transfType == 'trendPeaks') {
    cat(paste0('\nevaluating long term variations of the peaks'))
    TrendTh=try(tsEvaFindTrendThreshold(series, timeStamps, timeWindow),T)
    if(length(TrendTh)==0){
      TrendTh=NA
    }
    trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
    gevMaxima = 'annual'
    potEventsPerYear = epy
    minEventsPerYear = 1

  } else  if (transfType == 'trendChangeCIPercentile'){
    if (is.na(ciPercentile)){
      print('For trendCIPercentile transformation the label parameter cipercentile is mandatory')
    }
    cat('\n evaluating long term variations of extremes using the ', ciPercentile, 'th percentile and change point detection')
    trasfData = tsEvaTransformSeriesToStationaryTrendAndChangepts_ciPercentile(timeStamps, series, timeWindow,ciPercentile)
    gevMaxima = 'annual';
    potEventsPerYear = epy;
    minEventsPerYear = 0

  } else if (transfType == 'seasonalCIPercentile') {
    if (is.na(ciPercentile)) cat('For seasonalCIPercentile transformation the label parameter cipercentile is mandatory')
    cat(paste0('\nevaluating long term variations of extremes using the ', ciPercentile, 'th percentile\n'))
    trasfData = tsEvaTransformSeriesToStatSeasonal_ciPercentile( timeStamps, series, timeWindow, ciPercentile)
    gevMaxima = 'monthly'
    potEventsPerYear = 12
    minEventsPerYear = 6
  }
  if (args$potEventsPerYear != -1) potEventsPerYear = args$potEventsPerYear
  if (args$minEventsPerYear != -1) minEventsPerYear = args$minEventsPerYear

  if (dt<1) {
    pace=1/dt
    tsDaily=seq(1,length(trasfData$timeStamps),by=pace)
    trasfData$stdDevSeriesOr=trasfData$stdDevSeries
    trasfData$trendSeriesOr=trasfData$trendSeries
    trasfData$stdDevErrorOr=trasfData$stdDevError
    trasfData$stdDevSeries=trasfData$stdDevSeries[tsDaily]
    trasfData$trendSeries=trasfData$trendSeries[tsDaily]
    trasfData$stdDevError=trasfData$stdDevError[tsDaily]
  }

  ms = data.frame(trasfData$timeStamps, trasfData$stationarySeries)
  minPeakDistance = minPeakDistanceInDays/dt;

  #estimating the non stationary EVA parameters
  cat('\nExecuting stationary eva')
  pointData = tsEvaSampleData(ms, potEventsPerYear, minEventsPerYear, minPeakDistanceInDays,mode);
  evaAlphaCI = .68; # in a gaussian approximation alphaCI~68% corresponds to 1 sigma confidence
  eva = tsEVstatistics(pointData, evaAlphaCI, gevMaxima, gevType, evdType,shape_bnd);

  if (eva$isValid==FALSE) {
    print("problem in the computation of EVA statistics")
  }

  eva[[2]]$GPDstat$thresholdError <- pointData$POT$thresholdError

  # !!! Assuming a Gaussian approximation to compute the standard errors for
  # the GEV parameters
  if (eva[[2]][[1]]$method[1]!="No fit") {
    epsilonGevX <- eva[[2]][[1]]$parameters[3]
    errEpsilonX <- epsilonGevX - eva[[2]][[1]]$paramCIs[1,3]
    sigmaGevX <- eva[[2]][[1]]$parameters[2]
    errSigmaGevX <- sigmaGevX - eva[[2]][[1]]$paramCIs[1, 2]
    muGevX <- eva[[2]][[1]]$parameters[1]
    errMuGevX <- muGevX - eva[[2]][[1]]$paramCIs[1, 1]

    cat('\nTransforming to non stationary eva ...\n')
    epsilonGevNS = epsilonGevX;
    errEpsilonGevNS = errEpsilonX;
    sigmaGevNS = trasfData$stdDevSeries*sigmaGevX;

    #propagating the errors on stdDevSeries and sigmaGevX to sigmaGevNs.
    # err(sigmaNs) = sqrt{ [sigmaX*err(stdDev)]^2 + [stdDev*err(sigmaX)]^2 }
    # the error on sigmaGevNs is time dependant.
    errSigmaGevFit = trasfData$stdDevSeries*errSigmaGevX;
    errSigmaGevTransf = sigmaGevX*trasfData$stdDevError;
    errSigmaGevNS = (  errSigmaGevTransf^2   +  errSigmaGevFit^2  )^.5;
    muGevNS = trasfData$stdDevSeries*muGevX + trasfData$trendSeries;

    # propagating the errors on stdDevSeries, trendSeries and sigmaGevX to muGevNS.
    # err(muNs) = sqrt{ [muX*err(stdDev)]^2 + [stdDev*err(muX)]^2 + err(trend)^2 }
    # the error on muGevNS is time dependant.
    errMuGevFit = trasfData$stdDevSeries*errMuGevX;
    errMuGevTransf = (  (muGevX*trasfData$stdDevError)^2 + trasfData$trendError^2  )^.5;
    errMuGevNS = (  errMuGevTransf^2   +  errMuGevFit^2  )^.5;
    gevParams=c()
    gevParams$epsilon = epsilonGevNS;
    gevParams$sigma = sigmaGevNS;
    gevParams$mu = muGevNS;
    gevParams$annualMax=trasfData$nonStatSeries[pointData$annualMaxIndx]
    gevParams$monthlyMax=trasfData$nonStatSeries[pointData$monthlyMaxIndx]
    gevParams$annualMaxIndx=pointData$annualMaxIndx
    gevParams$monthlyMaxIndx=pointData$monthlyMaxIndx


    if(tolower(gevMaxima) == "annual") {
      gevParams$timeDelta <- 365.25
      gevParams$timeDeltaYears <- 1
    } else if(tolower(gevMaxima) == "monthly") {
      gevParams$timeDelta <- 365.25/12
      gevParams$timeDeltaYears <- 1/12
    }
    gevParamStdErr=c()
    gevParamStdErr$epsilonErr <- errEpsilonGevNS

    gevParamStdErr$sigmaErrFit <- errSigmaGevFit
    gevParamStdErr$sigmaErrTransf <- errSigmaGevTransf
    gevParamStdErr$sigmaErr <- errSigmaGevNS

    gevParamStdErr$muErrFit <- errMuGevFit
    gevParamStdErr$muErrTransf <- errMuGevTransf
    gevParamStdErr$muErr <- errMuGevNS

    gevObj=list()
    gevObj$method <- eva[[2]][[1]]$method
    gevObj$parameters <- gevParams
    gevObj$paramErr <- gevParamStdErr
    gevObj$stationaryParams <- eva[[2]][[1]]
    gevObj$objs$monthlyMaxIndexes <- pointData$monthlyMaxIndexes
  }else{

    epsilonGevX <- eva[[2]][[1]]$parameters[3]
    sigmaGevX <- eva[[2]][[1]]$parameters[2]
    muGevX <- eva[[2]][[1]]$parameters[1]
    cat('\nTransforming to non stationary eva ...\n')
    epsilonGevNS = epsilonGevX;
    sigmaGevNS = trasfData$stdDevSeries*sigmaGevX;
    muGevNS = trasfData$stdDevSeries*muGevX + trasfData$trendSeries;

    gevParams=c()
    gevParams$epsilon = epsilonGevNS;
    gevParams$sigma = sigmaGevNS;
    gevParams$mu = muGevNS;
    gevParams$annualMax=trasfData$nonStatSeries[pointData$annualMaxIndx]
    gevParams$monthlyMax=trasfData$nonStatSeries[pointData$monthlyMaxIndx]
    gevParams$annualMaxIndx=pointData$annualMaxIndx
    gevParams$monthlyMaxIndx=pointData$monthlyMaxIndx

    if(tolower(gevMaxima) == "annual") {
      gevParams$timeDelta <- 365.25
      gevParams$timeDeltaYears <- 1
    } else if(tolower(gevMaxima) == "monthly") {
      gevParams$timeDelta <- 365.25/12
      gevParams$timeDeltaYears <- 1/12
    }

    gevObj=list()
    gevObj$method = "No fit";
    gevObj$parameters = gevParams;
    gevObj$paramErr = NULL;
    gevObj$stationaryParams = NULL;
    gevObj$objs.monthlyMaxIndexes = NULL;
  }

  # estimating the non stationary GPD parameters
  # !!! Assuming a Gaussian approximation to compute the standard errors for
  # the GPD parameters
  if (eva[[2]][[2]]$method!="No fit") {
    epsilonPotX <- eva[[2]][[2]]$parameters[2]
    errEpsilonPotX <- epsilonPotX - eva[[2]][[2]]$paramCIs[1,1]
    sigmaPotX <- eva[[2]][[2]]$parameters[1]
    errSigmaPotX <- sigmaPotX - eva[[2]][[2]]$paramCIs[1, 2]
    thresholdPotX = eva[[2]][[2]]$parameters[3];
    errThresholdPotX = eva[[2]][[2]]$thresholdError;
    nPotPeaks = eva[[2]][[2]]$parameters[5];
    percentilePotX = eva[[2]][[2]]$parameters[6];

    dtPeaks = minPeakDistance;
    timeStamps=as.Date(timeStamps)
    dtPotX = as.numeric(timeStamps[length(timeStamps)] - timeStamps[1])/length(series)*dtPeaks;
    epsilonPotNS = epsilonPotX;
    errEpsilonPotNS = errEpsilonPotX;
    sigmaPotNS = sigmaPotX*trasfData$stdDevSeries;

    # propagating the errors on stdDevSeries and sigmaPotX to sigmaPotNs.
    # err(sigmaNs) = sqrt{ [sigmaX*err(stdDev)]^2 + [stdDev*err(sigmaX)]^2 }
    # the error on sigmaGevNs is time dependant.
    errSigmaPotFit = trasfData$stdDevSeries*errSigmaPotX;
    errSigmaPotTransf = sigmaPotX*trasfData$stdDevError;
    errSigmaPotNS = (  errSigmaPotTransf^2   +  errSigmaPotFit^2  )^.5;
    thresholdPotNS = thresholdPotX*trasfData$stdDevSeries + trasfData$trendSeries;
    # propagating the errors on stdDevSeries and trendSeries to thresholdPotNs.
    # err(thresholdPotNs) = sqrt{ [thresholdPotX*err(stdDev)]^2 + err(trend)^2 }
    # the error on thresholdPotNs is constant.
    thresholdErrFit = 0;

    thresholdErrTransf = ((trasfData$stdDevSeries*errThresholdPotX)^2 + (thresholdPotX*trasfData$stdDevError)^2  +  trasfData$trendError^2)^.5;
    thresholdErr = thresholdErrTransf;

    potParams=c()
    potParams$epsilon = epsilonPotNS;
    potParams$sigma = sigmaPotNS;
    potParams$threshold = thresholdPotNS;
    potParams$percentile = percentilePotX;
    potParams$timeDelta = dtPotX;
    potParams$timeDeltaYears = dtPotX/365.25;
    potParams$timeHorizonStart = min(trasfData$timeStamps);
    potParams$timeHorizonEnd = max(trasfData$timeStamps);
    potParams$peaks=trasfData$nonStatSeries[pointData$POT$ipeaks]
    potParams$peakID=pointData$POT$ipeaks
    potParams$peakST=pointData$POT$stpeaks
    potParams$peakEN=pointData$POT$endpeaks
    potParams$nPeaks = nPotPeaks;


    potParamStdErr=c()
    potParamStdErr$epsilonErr = errEpsilonPotNS;
    potParamStdErr$sigmaErrFit = errSigmaPotFit;
    potParamStdErr$sigmaErrTransf = errSigmaPotTransf;
    potParamStdErr$sigmaErr = errSigmaPotNS;
    potParamStdErr$thresholdErrFit = thresholdErrFit;
    potParamStdErr$thresholdErrTransf = thresholdErrTransf;
    potParamStdErr$thresholdErr = thresholdErr;

    potObj=list()
    potObj$method = eva[[2]][[2]]$method;
    potObj$parameters = potParams;
    potObj$paramErr = potParamStdErr;
    potObj$stationaryParams = eva[[2]][[2]];
    potObj$objs = NULL;
  }else{

    dtPeaks = minPeakDistance;
    timeStamps=as.Date(timeStamps)
    dtPotX = as.numeric(timeStamps[length(timeStamps)] - timeStamps[1])/length(series)*dtPeaks;
    thresholdPotX = pointData$POT$threshold
    thresholdPotNS = thresholdPotX*trasfData$stdDevSeries + trasfData$trendSeries;

    epsilonPotX <- pointData$POT$pars[2]
    sigmaPotX <- pointData$POT$pars[1]
    epsilonPotNS = epsilonPotX;
    sigmaPotNS = sigmaPotX*trasfData$stdDevSeries;
    thresholdPotNS = thresholdPotX*trasfData$stdDevSeries + trasfData$trendSeries;

    potParams=c()
    potParams$epsilon = epsilonPotNS;
    potParams$sigma = sigmaPotNS;
    potParams$threshold = thresholdPotNS;
    potParams$percentile = pointData$POT$percentile;
    potParams$timeDelta = dtPotX;
    potParams$timeDeltaYears = dtPotX/365.2425;
    potParams$timeHorizonStart = min(trasfData$timeStamps);
    potParams$timeHorizonEnd = max(trasfData$timeStamps);
    potParams$peaks=trasfData$nonStatSeries[pointData$POT$ipeaks]
    potParams$peakID=pointData$POT$ipeaks
    potParams$peakST=pointData$POT$stpeaks
    potParams$peakEN=pointData$POT$endpeaks
    potParams$nPeaks = length(pointData$POT$peaks);

    potObj=list()
    potObj$method = "No fit";
    potObj$parameters = potParams;
    potObj$paramErr = NULL;
    potObj$stationaryParams = NULL;
    potObj$objs = NULL;
  }

  # setting output objects
  nonStationaryEvaParams <- list(gevObj=gevObj, potObj=potObj)
  stationaryTransformData <- trasfData
  return(list(nonStationaryEvaParams=nonStationaryEvaParams,stationaryTransformData=stationaryTransformData))

}
