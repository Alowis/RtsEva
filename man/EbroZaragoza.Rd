% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{EbroZaragoza}
\alias{EbroZaragoza}
\title{Simulated river discharge of the Ebro river at Zaragoza}
\format{
Two column dataframe:
#' \itemize{
\item time (POSIXct timestamp with 6-hourly resolution)
\item Q (6-houry mean discharge in cubic meter per second)
}
}
\source{
\url{https://data.jrc.ec.europa.eu/dataset/a605a675-9444-4017-8b34-d66be5b18c95}
}
\usage{
data(EbroZaragoza)
}
\description{
A time series of simulated river discharge of the Ebro river at
Zaragoza (longitude = -0.825 \ latitude = 41.608) from
1951 to 2020. Time series extracted from the HERA dataset:
https://data.jrc.ec.europa.eu/dataset/a605a675-9444-4017-8b34-d66be5b18c95.
The Ebro is Spain's longest river, with low and high water levels
alternating throughout the year, influenced by winter snowmelt and summer
evaporation/human usage. The river is vital for agriculture.
}
\keyword{datasets}
