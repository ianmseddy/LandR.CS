globalVariables(c(
  ".", "age", "logAge", "pixelGroup", "..cohortDefinitionCols", "..addedColumns", "..neededCols", ":="
))

#'  Calculate climate effect
#'
#' Predict biomass change with climate variables
#'

#' @param cohortData The LandR cohortData object
#' @param pixelGroupMap the pixelGroupMap needed to match cohorts with raster values
#' @param cceArgs a list of datasets used by the climate function
#' @param year time of simulation - used to select from list of projected climate rasters
#' @param gmcsGrowthLimits lower and upper limits to the effect of climate on growth
#' @param gmcsMortLimits lower and upper limits to the effect of climate on mortality
#' @param gmcsMinAge minimum age for which to predict full effect of growth/mortality -
#'                   younger ages are weighted toward a null effect with decreasing age
#' @param cohortDefinitionCols cohortData columns that determine individual cohorts
#' @importFrom data.table setkey data.table
#' @importFrom reproducible Copy
#' @importFrom LandR asInteger
#' @importFrom terra compareGeom ncell
#' @importFrom stats na.omit predict median
#' @rdname calculateClimateEffect
#' @export
calculateClimateEffect <- function(cohortData, pixelGroupMap, cceArgs, year,
                                   gmcsGrowthLimits, gmcsMortLimits, gmcsMinAge,
                                   cohortDefinitionCols = c("age", "speciesCode", "pixelGroup")) {

   # extract relevant args

  #assume that climate normals are present in cceArgs as cceArts$historicalClimate$historical_<var>_normal
  #assume climate variables are denoted by cceArgs$climateVariables
  #
  anomalyVariables <- setdiff(names(climateVariables), "")
  allClimVar <- c(climateVariables, anomalyVariables)

  #get climate variables for pixelGroupMap

  if (length(anomalies) > 0) {
    #check if year is already subset - if not, subset to PSP period
    anomalyVar <-
    setnames(anomalyData, old = anomalies, new = names(anomalies))
    PSPmodelData <- anomalyData[PSPmodelData, on = c("OrigPlotID1")]

    #recalculate the anomaly(s) via subtraction
    for (i in 1:length(anomalies)) {
      #index because the name isn't preserved if you take the object itself
      temp <- anomalies[i]
      setnames(PSPmodelData, c(temp, names(temp)), c("var", "anom"))
      PSPmodelData[, anom := var - anom]
      newOrderCols <- names(PSPmodelData)[!names(PSPmodelData) %in% c("var", "anom") ]
      setcolorder(PSPmodelData, newOrderCols)
      setnames(PSPmodelData, c("var", "anom"), c(temp, names(temp)))
    }
  }





  ATA <- cceArgs$projectedClimateRasters[["ATA"]][[paste0("year", year)]]
  CMI <- cceArgs$projectedClimateRasters[["CMI"]][[paste0("year", year)]]
  CMInormal <- cceArgs$historicalClimateRasters[["CMI_normal"]]
  mcsModel <- cceArgs$mcsModel
  gcsModel <- cceArgs$gcsModel

  cohortData <- Copy(cohortData)
  neededCols <- c(cohortDefinitionCols, "B")
  neededCols <- neededCols[neededCols %in% colnames(cohortData)]
  climCohortData <- cohortData[, ..neededCols]

  if (!compareGeom(CMI, CMInormal, ATA, stopOnError = FALSE)) {
    stop("different number of pixels in the climate data. Please review how these are created")
  }


  CMIvals <- as.vector(CMI)
  CMInormalvals <- as.vector(CMInormal)
  ATAvals <- as.vector(ATA)
  pixels <- as.vector(pixelGroupMap)

  # Center observations on mean of original model data
  climateMatch <- data.table(
    "pixelGroup" = pixels,
    "CMI" = CMIvals,
    "ATA" = ATAvals,
    "CMInormal" = CMInormalvals
  )

  climateMatch <- climateMatch[!is.na(pixelGroup)]
  # Not all pixelGroups are in pixelGroupMap, because climCohortData is a subset
  # Take the median climate for each pixel group as some pixelgroups occur across multiple climate raster pixels
  climValues <- climateMatch[, .(
    "CMI" = median(CMI, na.rm = TRUE),
    "ATA" = median(ATA, na.rm = TRUE),
    "CMInormal" = median(CMInormal, na.rm = TRUE)
  ), by = "pixelGroup"]

  climCohortData[, logAge := log(age)]
  # set age = 0 to 1, to prevent -inf in prediction - this shouldn't affect predictions due to minimum age
  climCohortData[age == 0, logAge := 0]
  setkey(climCohortData, pixelGroup)
  setkey(climValues, pixelGroup)

  # Join cohort Data with climate data
  predData <- climCohortData[climValues]

  # remove NA values that exist only because of pixelGroupMap
  predData <- na.omit(predData)

  pixelGroupsPostSubset <- predData$pixelGroup
  agePostSubset <- predData$age
  speciesCodePostSubset <- predData$speciesCode

  modCohortDef <- FALSE

  # necessary for joining if cohortData has added columns
  if (length(cohortDefinitionCols[!cohortDefinitionCols %in% c("age", "pixelGroup", "speciesCode")]) > 0) {
    modCohortDef <- TRUE
    addedColumns <- cohortDefinitionCols[!cohortDefinitionCols %in% c("age", "pixelGroup", "speciesCode")]
    addedColumns <- predData[, ..addedColumns]
  }

  predData <- predData[, .(logAge, ATA, CMI, CMInormal)]

  # Create the 'reference climate' dataset to normalize the prediction
  refClim <- predData
  refClim$CMI <- refClim$CMInormal # replace CMI with the CMI normal for 1950-2010
  refClim$ATA <- 0 # the anomaly by definition has 0 as nromal

  refClim[, CMInormal := NULL] # or the mortality model will be upset
  predData[, CMInormal := NULL]

  # make growth prediction as ratio
  growthPred <- asInteger(predict(gcsModel, predData, level = 0, asList = TRUE, type = "response") /
                            predict(gcsModel, refClim, level = 0, asList = TRUE, type = "response") * 100)
  growthPred[growthPred < min(gmcsGrowthLimits)] <- min(gmcsGrowthLimits)
  growthPred[growthPred > max(gmcsGrowthLimits)] <- max(gmcsGrowthLimits)
  growthPred[is.na(growthPred)] <- max(gmcsGrowthLimits)

  # make mortality prediction

  mortPred <- asInteger(predict(
    object = mcsModel, parameter = "mu",
    newdata = predData, level = 0, asList = TRUE, type = "response"
  ) /
    predict(
      object = mcsModel, parameter = "mu", newdata = refClim,
      level = 0, asList = TRUE, type = "response"
    ) * 100)
  mortPred[mortPred < min(gmcsMortLimits)] <- min(gmcsMortLimits)
  mortPred[mortPred > max(gmcsMortLimits)] <- max(gmcsMortLimits)

  if (anyNA(c(mortPred, growthPred))) {
    mortPred[is.na(mortPred)] <- max(gmcsMortLimits)

    warning("NA in climate prediction. Likely integer overflow - setting to gmcsLimits")
  }

  # predict requires exact same columns in data.frame at the moment, hence this clumsy rebuilding
  climateEffect <- data.table(
    "pixelGroup" = pixelGroupsPostSubset,
    "speciesCode" = speciesCodePostSubset,
    "age" = agePostSubset,
    "growthPred" = growthPred,
    "mortPred" = mortPred
  )
  if (modCohortDef) {
    climateEffect <- cbind(climateEffect, addedColumns)
  }

  # restrict predictions to those above min stand age
  climateEffect[age < gmcsMinAge, growthPred := as.integer(100 + ((growthPred - 100) * (age / gmcsMinAge)))]
  climateEffect[age < gmcsMinAge, mortPred := as.integer(100 + ((mortPred - 100) * (age / gmcsMinAge)))]
  temp <- cohortData[, ..cohortDefinitionCols]
  climateEffect <- climateEffect[temp, on = cohortDefinitionCols]
  rm(temp)
  # this is to fix any pixelGroups that were dropped by the na.omit of climData due to NA climate values
  # which should be quite rare but persist with postProcess problems
  climateEffect[is.na(growthPred), c("growthPred", "mortPred") := .(100, 100)]

  # Because the params are numeric (e.g 66.667, the comparison forces the int to numeric)
  climateEffect[, c("growthPred", "mortPred") := .(asInteger(growthPred), asInteger(mortPred))]

  return(climateEffect)
}
