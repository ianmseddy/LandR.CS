test_that("getModelPred returns one prediction per pixelGroup/species", {
  cd <- data.table::data.table(
    pixelGroup   = c(1L, 1L, 2L),
    speciesCode  = factor(c("Pice_mar", "Pinu_ban", "Pice_mar"),
                          levels = c("Pice_mar", "Pinu_ban")),
    logAge       = c(log(60), log(60), log(40)),
    standBiomass = c(300, 300, 150),
    biomass      = c(100, 200, 150)
  )
  clim <- data.table::data.table(pixelGroup = c(1L, 2L), MAT = c(2, 3), CMI = c(1, 2))

  out <- getModelPred(cohortData = cd, climateCovariates = clim,
                      gmcsModel = testModel(), stat = "mortality")

  expect_s3_class(out, "data.table")
  expect_named(out, c("speciesCode", "pixelGroup", "mortality"), ignore.order = TRUE)
  expect_equal(nrow(out), 3L)
  expect_setequal(paste(out$pixelGroup, out$speciesCode), c("1 Pice_mar", "1 Pinu_ban", "2 Pice_mar"))
})

test_that("getModelPred predictions follow the fitted model", {
  cd <- data.table::data.table(
    pixelGroup   = c(1L, 2L),
    speciesCode  = factor(c("Pice_mar", "Pice_mar"), levels = c("Pice_mar", "Pinu_ban")),
    logAge       = c(2, 5),
    standBiomass = c(100, 200),
    biomass      = c(100, 200)
  )
  clim <- data.table::data.table(pixelGroup = c(1L, 2L), MAT = c(2, 3))

  ## intercept 1 + 2 * logAge, no biomass term => 5 and 11
  out <- getModelPred(cohortData = cd, climateCovariates = clim,
                      gmcsModel = testModel(intercept = 1, coefLogAge = 2),
                      stat = "mortality")
  data.table::setkey(out, pixelGroup)
  expect_equal(out$mortality, c(5, 11))
})

test_that("growth predictions are exponentiated but mortality predictions are not", {
  cd <- data.table::data.table(
    pixelGroup   = 1L,
    speciesCode  = factor("Pice_mar", levels = c("Pice_mar", "Pinu_ban")),
    logAge       = 3,
    standBiomass = 100,
    biomass      = 100
  )
  clim <- data.table::data.table(pixelGroup = 1L, MAT = 2)
  mod <- testModel()

  growth <- getModelPred(cd, clim, gmcsModel = mod, stat = "growth")
  mort   <- getModelPred(cd, clim, gmcsModel = mod, stat = "mortality")

  expect_equal(growth$growth, exp(mort$mortality))
})

test_that("getModelPred averages across folds", {
  cd <- data.table::data.table(
    pixelGroup   = 1L,
    speciesCode  = factor("Pice_mar", levels = c("Pice_mar", "Pinu_ban")),
    logAge       = 2,
    standBiomass = 100,
    biomass      = 100
  )
  clim <- data.table::data.table(pixelGroup = 1L, MAT = 2)

  ## Fold1 predicts logAge (= 2); Fold2 predicts 10 + logAge (= 12). Mean = 7.
  twoFolds <- list(Fold1 = testFold(intercept = 0), Fold2 = testFold(intercept = 10))
  out <- getModelPred(cd, clim, gmcsModel = twoFolds, stat = "mortality")

  expect_equal(out$mortality, 7)
})

test_that("species in the fitted model but absent from cohortData are added as zeroes", {
  cd <- data.table::data.table(
    pixelGroup   = 1L,
    speciesCode  = factor("Pice_mar", levels = c("Pice_mar", "Pinu_ban")),
    logAge       = 2,
    standBiomass = 100,
    biomass      = 100
  )
  clim <- data.table::data.table(pixelGroup = 1L, MAT = 2)

  ## The fitted model knows about a species this landscape does not contain.
  mod <- list(Fold1 = testFold(extraCols = "sppBetu_pap"),
              Fold2 = testFold(extraCols = "sppBetu_pap"))

  out <- getModelPred(cd, clim, gmcsModel = mod, stat = "mortality")

  ## The absent species is filled with 0 (not present), so it must not appear in
  ## the output, and the present species must still be predicted.
  expect_equal(nrow(out), 1L)
  expect_equal(as.character(out$speciesCode), "Pice_mar")
})

test_that("getModelPred rejects an unknown stat", {
  cd <- data.table::data.table(pixelGroup = 1L, speciesCode = factor("Pice_mar", levels = c("Pice_mar", "Pinu_ban")),
                               logAge = 2, standBiomass = 100, biomass = 100)
  clim <- data.table::data.table(pixelGroup = 1L, MAT = 2)

  expect_error(getModelPred(cd, clim, gmcsModel = testModel(), stat = "biomass"))
})
