
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RtsEva

<!-- badges: start -->

[![R-CMD-check](https://github.com/Alowis/RtsEva/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Alowis/RtsEva/actions/workflows/R-CMD-check.yaml)
[![PRsWelcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](https://makeapullrequest.com/)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/RtsEva)](https://cran.r-project.org/package=RtsEva)
![](https://cranlogs.r-pkg.org/badges/grand-total/RtsEva)

<!-- badges: end -->

<img src="man/figures/RtsEVA.png" align="right" height="200"/>

This package is an adaptation of the Matalb tsEVA toolbox developed by
Lorenzo Mentaschi availaible here: <https://github.com/menta78/tsEva>

It contains an implementation of the Transformed-Stationary (TS)
methodology for non-stationary EVA as described in Mentaschi et
al. (2016). In synthesis this approach consists in (i) transforming a
non-stationary time series into a stationary one to which the stationary
EVA theory can be applied; and (ii) reverse-transforming the result into
a non-stationary extreme value distribution.

## References

[Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L.,
Besio, G., and Alfieri, L.: The transformed-stationary approach: a
generic and simplified methodology for non-stationary extreme value
analysis, Hydrol. Earth Syst. Sci., 20,3527-3547,
doi:10.5194/hess-20-3527-2016,
2016](http://www.hydrol-earth-syst-sci.net/20/3527/2016/)

## Installation

You can install the released version of mobirep from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RtsEva")
```
Alternatively, you can install the development version of RtsEva from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Alowis/RtsEva")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(RtsEva)
# Load a time series
timeAndSeries <- ArdecheStMartin
# go from six-hourly values to daily max
timeAndSeries <- max_daily_value(timeAndSeries)

# set a temporal window for the computation of running statistics
timeWindow <- 30*365 # 30 years

# Run the non-stationnary EVA
result <- TsEvaNs(timeAndSeries, timeWindow,
transfType = 'trendPeaks',tail = 'high')
```

After fitting the non-stationnay EVA, the package offers functions to
visualize the plots

## More

https://alowis.quarto.pub/rtseva-package-demo/

## Contact

For any questions or inquiries, please contact the package maintainer at
<alois.tilloy@ec.europa.eu>
