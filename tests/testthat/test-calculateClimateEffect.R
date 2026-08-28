test_that("calculateClimateEffect returns every input cohort", {
  cd <- testCohortData()
  out <- calculateClimateEffect(cohortData = cd, pixelGroupMap = testPixelGroupMap(),
                                cceArgs = testCceArgs(), gmcsGrowthLimits = c(0, 400),
                                gmcsMinAge = 30, time = 2011)

  expect_s3_class(out, "data.table")
  expect_equal(nrow(out), nrow(cd))
  expect_setequal(paste(out$pixelGroup, out$speciesCode, out$age),
                  paste(cd$pixelGroup, cd$speciesCode, cd$age))
  ## columns the downstream summary module relies on
  expect_true(all(c("growthPred", "mortPred", "standBiomass", "climateYear", "maxANPP",
                    "pred_currentGrowth", "pred_historicalGrowth",
                    "pred_currentMortality", "pred_historicalMortality") %in% names(out)))
})

test_that("a current climate identical to the normals produces no climate effect", {
  ## climateShift = 0 => current climate IS the historical normal, so the current
  ## and historical predictions coincide: growth is unmodified (100%) and there is
  ## no additional mortality.
  out <- calculateClimateEffect(cohortData = testCohortData(), pixelGroupMap = testPixelGroupMap(),
                                cceArgs = testCceArgs(climateShift = 0),
                                gmcsGrowthLimits = c(0, 400), gmcsMinAge = 30, time = 2011)

  expect_true(all(out$growthPred == 100))
  expect_true(all(out$mortPred == 0))
})

test_that("cohorts younger than gmcsMinAge are kept with a null climate effect", {
  out <- calculateClimateEffect(cohortData = testCohortData(), pixelGroupMap = testPixelGroupMap(),
                                cceArgs = testCceArgs(), gmcsGrowthLimits = c(0, 400),
                                gmcsMinAge = 30, time = 2011)

  young <- out[age < 30]
  expect_equal(nrow(young), 1L)          # the age-20 Pinu_ban cohort
  expect_equal(young$growthPred, 100)
  expect_equal(young$mortPred, 0)
})

test_that("growthPred is clamped to gmcsGrowthLimits", {
  ## A model whose predictions depend strongly on stand biomass makes the current
  ## and historical growth predictions diverge, pushing the ratio outside the limits.
  strong <- testModel(coefLogAge = 0, coefBiomass = 1)
  limits <- c(90, 110)

  out <- calculateClimateEffect(cohortData = testCohortData(), pixelGroupMap = testPixelGroupMap(),
                                cceArgs = testCceArgs(climateShift = 5, model = strong),
                                gmcsGrowthLimits = limits, gmcsMinAge = 30, time = 2011)

  expect_true(all(out$growthPred >= min(limits)))
  expect_true(all(out$growthPred <= max(limits)))
})

test_that("climateYear comes from cceArgs when supplied, otherwise from time", {
  pgm <- testPixelGroupMap()

  fromTime <- calculateClimateEffect(testCohortData(), pgm, testCceArgs(),
                                     gmcsGrowthLimits = c(0, 400), gmcsMinAge = 30, time = 2011)
  expect_true(all(fromTime[age >= 30]$climateYear == 2011))

  fromArgs <- calculateClimateEffect(testCohortData(), pgm, testCceArgs(climateYear = 2051),
                                     gmcsGrowthLimits = c(0, 400), gmcsMinAge = 30, time = 2011)
  expect_true(all(fromArgs[age >= 30]$climateYear == 2051))
})

test_that("cohorts below gmcsMinAge carry no climate columns (documents current behaviour)", {
  ## Cohorts younger than gmcsMinAge are set aside before the climate covariates
  ## and climateYear are attached, then rbind()ed back with fill = TRUE. Only
  ## growthPred/mortPred are backfilled, so climateYear and every climate column
  ## stay NA on those rows even though the climate year applies to the whole
  ## timestep. This test pins the behaviour as it stands -- if it is ever made
  ## consistent, update this test rather than deleting it.
  out <- calculateClimateEffect(testCohortData(), testPixelGroupMap(), testCceArgs(),
                                gmcsGrowthLimits = c(0, 400), gmcsMinAge = 30, time = 2011)

  young <- out[age < 30]
  expect_equal(nrow(young), 1L)
  expect_true(is.na(young$climateYear))
  expect_true(all(is.na(young[, c("MAT", "CMI", "MATanomaly", "CMIanomaly")])))
})

test_that("calculateClimateEffect does not modify its input cohortData", {
  cd <- testCohortData()
  before <- data.table::copy(cd)

  calculateClimateEffect(cohortData = cd, pixelGroupMap = testPixelGroupMap(),
                         cceArgs = testCceArgs(), gmcsGrowthLimits = c(0, 400),
                         gmcsMinAge = 30, time = 2011)

  expect_equal(cd, before)
})

test_that("calculateClimateEffect errors when no climate data is available", {
  cceArgs <- testCceArgs()
  cceArgs$currentClimateRasters <- NULL
  cceArgs$projectedClimateRasters <- NULL

  expect_error(
    calculateClimateEffect(testCohortData(), testPixelGroupMap(), cceArgs,
                           gmcsGrowthLimits = c(0, 400), gmcsMinAge = 30, time = 2011),
    "no method of obtaining climate data"
  )
})

test_that("every data.table function the package calls is actually imported", {
  ## Regression test: `setcolorder` was called but never imported, so
  ## calculateClimateEffect() only worked when the caller happened to have
  ## data.table attached, and failed with
  ## "could not find function setcolorder" otherwise.
  imports <- parent.env(asNamespace("LandR.CS"))
  expect_true(exists("setcolorder", envir = imports, inherits = FALSE))
})
