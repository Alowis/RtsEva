
All notable changes to this project will be documented in this file.

## [RtsEva 1.2.0] - 2026-09-22

### Changed

* 'tsGetPOT' now selects the optimal POT threshold using a continuous
penalty on the deviance (based on the deficit from 'minEventsPerYear') together
with a shape-parameter boundary/normalized-distance fallback score, replacing
the previous AIC + skip/penalty scheme. It now fits the GPD with 'L-BFGS-B'
constrained by 'shape_bnd' (std.err.type = "observed") and returns the full fit
object in 'pars'.
* 'tsGetPOT' uses a finer grid of candidate percentiles above the 95th
percentile and cleans NA/Inf values before peak detection.
* 'shape_bnd' is now a user-facing argument of 'tsGetPOT', 'tsEvaSampleData'
and 'TsEvaNs'. In 'TsEvaNs' it defaults automatically to c(-0.5, 1) for the
high tail and c(-1, 0) for the low tail when left as NA.
* 'tsEVstatistics' now reuses the GPD fit computed during POT threshold
selection (stored in 'pointData$POT$pars') instead of re-fitting the GPD.
* 'TsEvaNs' gained additional robustness in the 'trendPeaks' case when the
automatic trend threshold cannot be estimated, and guards the low-flow
transformation against a NULL 'trans'.
* 'tsEvaFindTrendThreshold' keeps the stability, negative-flow and percentile
vectors aligned and guards the breaking-point detection against short series.

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
