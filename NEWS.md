version 2.0.0.9007
- makeClimateVariablesForGMCS now available to prepare climate variables for simulations and model fitting
version 2.0.0.9005
- all model columns and predictions are now returned, includign the historical and current growth and mortality and all climate covariates
- The function determines which climate data to use based on the contents of cceArgs, using the following precedence order to maintain backward compatibility with earlier interfaces:
  currentClimateRasters (highest precedence) - if cceArgs$currentClimateRasters is supplied, these rasters are used directly and no other climate arguments are evaluated.
  climateYear - if currentClimateRasters is not supplied and cceArgs$climateYear is provided, the function retrieves the corresponding year from projectedClimateRasters.
  Time - if neither currentClimateRasters nor climateYear are supplied, the function attempts to retrieve the raster matching Time from projectedClimateRasters
  Failure condition - if none of the above are available (or if projectedClimateRasters is missing when required), the function stops with an error.
 


version 2.0.0

# this represents a major overhaul
- The main change is to accomodate the new class and structure of statistical model from the latest gmcsDataPrep module. 
 Specifically, `gcsModel` and `gmsModel` objects in cceArgs are assumed to be xgboost models, with tree
 species as a covariate, as well as arbitrary climate covariates (including their anomalies) and stand biomass. 
other changes include the following:
- functions supporting gamlss predictions have been dropped as this was deprecated in favor of xgboost
- two options have been added to the package supporting output of predictions, `debug` and `logPath`
  If `debug` is TRUE, a `data.table` containing covariates and predictions will be written to the `logPath` directory.
  This will have a simple name of `gmcsPred` and the year the prediction is made. Speaking of year...
- `time` was added as an argument to `calculateClimateEffect`for the purpose of writing output files during debug mode.
- mortality is no longer treated as a multiplier; instead it will have an additive effect on cohortData. If multiple cohorts
are present, the predicted species mortality will be divided between them proportional to their biomass. 
- Due to the above, gmcsMortLimits was removed as a parameter.  
- The package can optionally handle climateYear, a numeric vector denoting a year within a climate stack, e.g. 2020, 
- climateYear could potentially be included in cceArgs when simulations require resampling of historical or projected data 
(e.g. spin ups, NRV), and can be supplied by running the module `PredictiveEcology/climateYear` 


Known issues: https://github.com/ianmseddy/LandR.CS/issues

version 0.0.1
=============

## bug fixes

* `LandR.CS::own()` no longer uses `grep(..., sys.calls())` which takes many many minutes if `...` is a simList. Use `reproducible::.grepSysCalls` which is for this purpose.


