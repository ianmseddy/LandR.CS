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
#' @importFrom data.table `:=` as.data.table data.table rbindlist set setkey setnames 
#' @importFrom LandR asInteger
#' @importFrom terra ncell rast
#' @importFrom stats na.omit predict median
#' @rdname calculateClimateEffect
#' @export
calculateClimateEffect <- function(cohortData, pixelGroupMap, cceArgs,
                                   gmcsGrowthLimits, gmcsMortLimits, gmcsMinAge,
                                   cohortDefinitionCols = c("age", "speciesCode", "pixelGroup")) {
  
  originalCD <- cohortData
  cohortData <- copy(cohortData)
   # extract relevant args
  
  #assume that climate normals are present in cceArgs as cceArts$historicalClimate$historical_<var>_normal
  #assume climate variables are denoted by cceArgs$climateVariables
  climateVariables <- cceArgs$climateVariablesForGMCS
  climateRasters <- cceArgs$projectedClimateRasters
   rasterToGet <- paste0("year", cceArgs$climateYear)
  #make composite raster
  climateRasters <- lapply(climateRasters[climateVariables], "[[", rasterToGet) |>
    rast()
  
  # Get the variables which require anomalies
  # Note if anomaly variables become available directly, they should be 'unnamed'
  # to skip the steps that manually calculate them
  anomalyVariables <- setdiff(names(climateVariables), "")
  if (length(anomalyVariables)) {
    anomalyVariables <- setdiff(names(climateVariables), "")
    varsWithAnomalies <- climateVariables[anomalyVariables]
    normalNames <- paste0(varsWithAnomalies, "_normal")
    climateNormal <- cceArgs$historicalClimateRasters[normalNames] |>
      rast()
    anomalyRasters <- terra::subset(climateRasters, varsWithAnomalies)
    names(climateNormal) <- names(anomalyRasters)
    #TODO: think of ways this might get screwed up - prevent it 
    anomalyOut <- lapply(varsWithAnomalies, 
                         function(anomaly, normals = climateNormal, current = anomalyRasters) {
                           anomalyRas <- current[[anomaly]] - normals[[anomaly]]
                         }) |>
      rast()
    names(anomalyOut) <- anomalyVariables
    allClimateRas <- c(climateRasters, anomalyOut)
  } else { 
    allClimateRas <- climateRasters 
  }
  
  climateCovariates <- as.data.table(allClimateRas, cells = TRUE)
  pixelGroupDT <- data.table(pixelIndex = 1:ncell(pixelGroupMap), pixelGroup = as.vector(pixelGroupMap))
  climateCovariates <- climateCovariates[pixelGroupDT, on = c("cell" = "pixelIndex")]
  #calculate the mean climate variable for each pixelGroup, as some groups overlap multiple pixels
  climateCovariates <- climateCovariates[!is.na(pixelGroup), lapply(.SD, mean), 
                    .SDcols = names(allClimateRas), .(pixelGroup)]
  
  historicalClimate <- as.data.table(climateNormal, cells = TRUE)
  historicalClimate <- historicalClimate[pixelGroupDT, on = c("cell" = "pixelIndex")]
  historicalClimate <- historicalClimate[!is.na(pixelGroup), lapply(.SD, mean), 
                                         .SDcol = climateVariables, .(pixelGroup)]
  #assume that anomalies are 0 (ie they are the normal)
  if (length(anomalyVariables)){
    set(historicalClimate, NULL, anomalyVariables, 0)
  }
  #TODO: 
  #dataprep should use a function to at least confirm colnames and units correct
  # standBiomass + biomass in g/m2, logAge = logged standAge (ie identical by cohort in PG),
 
  #subset cd for necessary columns
  cohortData <- cohortData[, .(speciesCode, B, mortality, age, pixelGroup)]
  
  #calculate biomass-weighted age, then log it
  cohortData[, standBiomass := sum(B), .(pixelGroup)]
  cohortData[, weightedAge := sum(B * age)/standBiomass, .(pixelGroup)]
  cohortData[, logAge := log(weightedAge)]
  #we can predict species but not individual cohorts (this is absent from PSP, as we do not have their ages)
  #in the event there are other columns in cohortData - subset like so
  cd <- cohortData[, .(biomass = sum(B)), .(pixelGroup, speciesCode, logAge, standBiomass)]

  #Make predictions
  predictedGrowth <- lapply(list(climateCovariates, historicalClimate), 
                            getModelPred, cohortData = cd, 
                            gmcsModel = cceArgs$gcsModel, 
                            stat = "growth")
  #take the 'historical' prediction for calculation of modifier
  subsetDT <- predictedGrowth[[2]][, .(pixelGroup, speciesCode, growth)]
  setnames(subsetDT, "growth", "historicalGrowth")
  #take means before join to avoid having to track fold

  predictedGrowth <- predictedGrowth[[1]][subsetDT, on = c("pixelGroup", "speciesCode")]
  predictedGrowth[, growthPred := asInteger(growth/historicalGrowth * 100)]
  cohortData <- cohortData[predictedGrowth, on = c("pixelGroup", "speciesCode")]
  # TODO: Reduce duplication and wrap these 2 in a function
  # this final wrapper is hurting my brain right now
  # getModelPred is needlessly building the model matrix 4 times (historical/current/growth/mort)
  # the outer wrapper should this?
  predictedMortality <- lapply(list(climateCovariates, historicalClimate), 
                            getModelPred, cohortData = cd, 
                            gmcsModel = cceArgs$mcsModel, 
                            stat = "mortality")
  #take the 'historical' prediction for calculation of modifier
  subsetDT <- predictedMortality[[2]][, .(pixelGroup, speciesCode, mortality)]
  setnames(subsetDT, "mortality", "historicalMortality")
  #take means before join to avoid having to track fold
  
  predictedMortality <- predictedMortality[[1]][subsetDT, on = c("pixelGroup", "speciesCode")]
  predictedMortality[, mortPred := asInteger(mortality/historicalMortality * 100)]
  cohortData <- cohortData[predictedMortality, on = c("pixelGroup", "speciesCode")]
  #could do this in getModelPred but I prefer to see the preds prior to ceiling/floor application
  cohortData[, mortPred := pmax(min(gmcsMortLimits), pmin(mortPred, max(gmcsMortLimits)))]
  cohortData[, growthPred := pmax(min(gmcsGrowthLimits), pmin(growthPred, max(gmcsGrowthLimits)))]
  
  #TODO: discuss how to deal with this - does it affect b-weighted standAge?
  cohortData[age < gmcsMinAge, .(growthPred, mortPred) := c(100, 100)]
  return(cohortData)
}

#' getModelPred
#'
#' Join cohortData with either the growth or mortality xgboost models
#' This is a utility function to avoid repetition as 2 predictions are needed
#' for each model (under the current climate, and historical) to achieve the 
#' climate modifier

#' @param cohortData The LandR cohortData object
#' @param climateCovariates data.table of mean climate covariates by pixelGroup 
#' @param stat character vector of either "growth" or "mortality"
#' this is used to name the pred column and exponentiate predicted growth if applicable
#' @importFrom data.table setkey data.table
#' @importFrom LandR asInteger
#' @importFrom terra ncell
#' @importFrom data.table as.data.table copy melt.data.table rbindlist setnames `.` `.SD`
#' @importFrom stats model.matrix
#' @keywords internal
#' @rdname getModelPred
getModelPred <- function(cohortData, climateCovariates, gmcsModel, stat = "growth") {
  
  stopifnot(stat %in% c("growth", "mortality"))
  
  predData <- cohortData[climateCovariates, on = "pixelGroup"]
  setnames(predData, "speciesCode", "spp") #to match fitting
  #take the character vector and make into dummy variables like so
  #give each species "spp" at the start to match the fitted data (ideally this is not hardcoded)
  predData <- model.matrix(~ . + 0, data = predData)
  predData <- as.data.table(predData)
  
  #if tree species are absent, they must be added for the predict
  fittedCols <- names(gmcsModel$Fold1$valData)
  missing <- setdiff(fittedCols, c(names(predData), "obs", "pred", "resid"))
  if (length(missing)){
    set(predData, j = missing, value =  0)
  }
  
  gmPreds <- lapply(gmcsModel, function(Fold, data = predData){
    out <- copy(data)
    out[, predCol := predict(object = Fold$mod, newdata = data)]
    return(out)
  }) |>
    rbindlist()
  #now summarize 
  #probably you will have to do this 4 times so turn it into a function
  out <- melt.data.table(gmPreds, id.vars = c("pixelGroup", "predCol"), 
                         measure.vars = grep("spp", colnames(gmPreds), value = TRUE),
                         variable.name = "speciesCode", value.name = "presence")
  out <- out[presence == 1,] #keep only presences
  out[, speciesCode := gsub(pattern = "spp", replacement = "", speciesCode)]
  out[, speciesCode := factor(speciesCode, levels = levels(cohortData$speciesCode))]
  
  #take the mean of the folds
  out <- out[, .(predCol = mean(predCol)), .(speciesCode, pixelGroup)]
  
  if (stat == "growth") {
    out[, predCol := exp(predCol)]
  }
  #rename column to either mortality or growth
  setnames(out, "predCol", stat)
  
  return(out)
}