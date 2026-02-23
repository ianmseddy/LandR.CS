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


