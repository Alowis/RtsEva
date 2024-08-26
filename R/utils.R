#Utilities used in other RTSEVA functions

#' Fill missing values in a time series using a moving average approach.
#'
#' This function takes a vector of timestamps and a corresponding series with missing values,
#' and fills the missing values by taking the average of the surrounding values.
#'
#' @param timeStamps A vector of timestamps.
#' @param series A vector representing the time series with missing values.
#'
#' @return A vector with missing values filled using a moving average approach.
#'
#' @examples
#' timeStamps <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' series <- c(1, 2, NA, 4, 5, NA, 7, 8, NA, 10)
#' filledSeries <- tsEvaFillSeries(timeStamps, series)
#' filledSeries
#'
#' @export
tsEvaFillSeries <- function(timeStamps, series) {
  indxs <- which(is.na(series))
  for (ix in indxs){
    id1=max(1,ix-4)
    id2=min(length(series),ix+4)
    instant=series[c(id1:id2)]
    filler=mean(instant,na.rm=T)
    series[ix]=filler
  }
  return(series)
}

#' Parse named arguments and assign values to a predefined argument structure.
#'
#' This function takes a list of named arguments and assigns their values to a predefined argument structure.
#' The argument structure is a list with named elements representing the available arguments.
#' If an argument is present in the list of named arguments, its value is assigned to the corresponding element in the argument structure.
#' If an argument is not present, its value in the argument structure remains unchanged.
#'
#' @param args A list of named arguments.
#' @param argStruct A list representing the argument structure with named elements.
#'
#' @return A modified argument structure with values assigned from the list of named arguments.
#'
#' @examples
#' args <- list(arg1 = 10, arg2 = "tanargue")
#' argStruct <- list(arg1 = 0, arg2 = "", arg3 = TRUE)
#' modifiedArgStruct <- tsEasyParseNamedArgs(args, argStruct)
#' modifiedArgStruct
#'
#' @export
tsEasyParseNamedArgs <- function(args, argStruct) {
  avlArgs <- names(argStruct)
  for (ia in 1:length(avlArgs)) {
    argName <- avlArgs[ia]
    argIndx <- which(names(args) == argName)
    if (length(argIndx)>0) {
      val <- args[[argIndx]]
      argStruct[[ia]] <- val
    }
  }
  return(argStruct)
}

#' Check if all years in a time series are present
#'
#' This function checks if all years specified in a given time series are present.
#'
#' @param timeseries A time series object.
#' @param yro A vector specifying the start and end years.
#'
#' @return A logical value indicating whether all years in the time series are present.
#'
#' @examples
# Create a time series object
#' ts_data <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2004-12-31"), by = "year")

# Check if all years in the time series are present
#' check_timeseries(ts_data, c(2000, 2004))
#' # Output: TRUE
#'
#' check_timeseries(ts_data, c(2000, 2005))
#' # Output: FALSE
#'
#' @importFrom lubridate year
#' @export
check_timeserie=function(timeseries,yro){
  ts_years <- as.integer((lubridate::year(timeseries)))
  year_check <- yro %in% ts_years
  runs <- rle(year_check)
  rf=which(runs$values==FALSE)
  if (any(runs$lengths[rf] >= 4)) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}

#' Max Daily Value Function
#'
#' This function converts a 6-hourly time series to a daily time series and
#' calculates the maximum value for each day.
#'
#' @param timeseries A time series with a 6-hourly resolution.
#' @return A data frame containing the daily maximum values.
#' @export
#' @importFrom xts apply.daily
#' @examples
#' # Example usage:
#'timeseries <- ArdecheStMartin
#' max_daily_value(timeseries)
max_daily_value <- function(timeseries) {
  # Convert the 6-hourly time series to daily time series
  daily_timeseries <- apply.daily(timeseries, max)
  tday=unique(as.Date(as.character(timeseries[,1])))
  return(data.frame(date=tday,Qmd=daily_timeseries))
}
