## Shared fixtures for LandR.CS tests.
##
## The real growth/mortality models are xgboost objects produced by the
## `gmcsDataPrep` module. The package never touches xgboost directly -- it only
## calls the `predict()` generic and reads `gmcsModel$Fold1$valData` for its
## column names. So a plain `lm` stands in perfectly, which keeps the tests fast
## and free of an xgboost dependency.

## A 4x4 raster with the geometry every fixture shares.
testRaster <- function(vals) {
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
  terra::values(r) <- vals
  r
}

## pixelGroupMap: 6 cells in group 1, 6 in group 2, 4 unassigned (NA).
testPixelGroupMap <- function() {
  testRaster(c(rep(1L, 6), rep(2L, 6), rep(NA_integer_, 4)))
}

## One "fold" of a fitted model. The training data is an exact linear function of
## its predictors, so the fit is perfect and predictions are deterministic:
##   pred = intercept + coefLogAge * logAge + coefBiomass * standBiomass
testFold <- function(intercept = 0, coefLogAge = 1, coefBiomass = 0, extraCols = character()) {
  trn <- data.frame(logAge = c(1, 2, 3, 4, 5, 6),
                    standBiomass = c(10, 35, 20, 60, 45, 80))
  trn$y <- intercept + coefLogAge * trn$logAge + coefBiomass * trn$standBiomass
  valData <- data.frame(logAge = 0, standBiomass = 0, obs = 0, pred = 0, resid = 0)
  for (nm in extraCols) valData[[nm]] <- 0
  list(mod = stats::lm(y ~ logAge + standBiomass, data = trn), valData = valData)
}

testModel <- function(...) list(Fold1 = testFold(...), Fold2 = testFold(...))

## cohortData with two pixelGroups and two species. The age-20 cohort sits below
## the default gmcsMinAge used in the tests, which exercises the "too young" path.
testCohortData <- function() {
  data.table::data.table(
    pixelGroup  = c(1L, 1L, 2L, 2L),
    speciesCode = factor(c("Pice_mar", "Pinu_ban", "Pice_mar", "Pinu_ban")),
    age         = c(50L, 80L, 60L, 20L),
    B           = c(100, 200, 150, 50),
    mortality   = c(1, 2, 3, 4),
    maxANPP     = c(5, 6, 7, 8)
  )
}

## cceArgs for the supported code path: every climate variable is *named*, so each
## one gets an anomaly layer computed against its historical normal.
##
## `climateShift` offsets the current climate away from the normals. At the
## default of 0 the current climate IS the normal, so current and historical
## predictions coincide -- the "no climate effect" case.
testCceArgs <- function(climateShift = 0, model = testModel(), ...) {
  normals <- list(MAT_normal = testRaster(rep(2, 16)),
                  CMI_normal = testRaster(rep(1, 16)))
  current <- c(testRaster(rep(2 + climateShift, 16)),
               testRaster(rep(1 + climateShift, 16)))
  names(current) <- c("MAT", "CMI")

  c(list(climateVariablesForGMCS = c(MATanomaly = "MAT", CMIanomaly = "CMI"),
         currentClimateRasters = current,
         historicalClimateRasters = normals,
         gcsModel = model,
         mcsModel = model), list(...))
}
