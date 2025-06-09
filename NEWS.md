# News

All notable changes to this project will be documented in this file.

## [RtsEva 1.1.0] - 2025-06-09

### Added

* 'tsEvaTransformSeriesToStationaryMMXTrend()' added. It computes the trend of
monthly maxima and adds to the suite of trend computation functions.


### Changed

* New rules for the selection of the optimal GPD fit in 'tsGetPOT'. The fit is 
now constrained based on the shape parameter value (need to be between two bounds)
and the AIC. 
* TrendTH can be specified outside the 'TsEvaNs' function
* Shape bounds updated in 'TsEvaNs'
* Updated handling of trendPeaks cases in 'TsEvaNs': increase robustness with 
iterative approach in cases where the trend is completely stable. 
* 'check_timeseries' now accepts timeseries where a maximum of two years are missing


### Fixed

* Correction of a mistake in the output writing of 'tsEvaComputeReturnLevelsGEV'. 
The output matrix was not initiated properly. The new matrix is a transposition 
the old one.


## [1.0.0] - 2024-06-24
