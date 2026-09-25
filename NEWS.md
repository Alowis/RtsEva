
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
* 'tsEvaNanRunningMean' was reimplemented with cumulative sums for a large
speed-up on long time series; results are unchanged. Package byte-compilation
was enabled ('ByteCompile: true').

### Fixed

* 'tsEvaNanRunningVariance' returned a slightly incorrect running variance: in
the previous incremental implementation the count of valid points could drift
out of sync with the summed squared values (the removal used index
'minindx - 1' while the addition used 'maxindx + 1'). It has been reimplemented
with cumulative sums to compute the correct centered-window mean of squares.
The correction to the variance is small: on the bundled ArdecheStMartin series,
the relative change is of the order of 1-2% (median ~1.7% at a one-year window),
so fitted GEV/GPD standard-deviation-dependent parameters shift only slightly.
* 'tsEvaNanRunningStatistics' contained more serious errors: the same
valid-point count desynchronisation and, in addition, each incoming value was
centered by the running mean at an incorrect index. This produced third and
fourth running moments that were substantially wrong (median relative errors of
several hundred percent on the bundled ArdecheStMartin series). It has been
reimplemented with cumulative sums to compute the correct centered-window
moments (validated to machine precision against a direct definition). Diagnostic
quantities derived from these moments therefore change from previously incorrect
values to correct ones.
* As a consequence of the running-variance fix,
'tsEvaTransformSeriesToStatSeasonal_ciPercentile' no longer returns an all-NA
'trendSeries' on short series. Previously the incorrect running variance could
yield negative values whose square root produced NaNs that propagated through
the seasonal standard-deviation estimation; this is now resolved.

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
