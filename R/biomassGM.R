#' Predict biomass change with climate variables
#'

#' @param cohortData The LandR cohortData object
#' @param pixelGroupMap the pixelGroupMap needed to match cohorts with raster values
#' @param cceArgs a list of additional data and arguments used by the climate function
#' @param gmcsGrowthLimits lower and upper limits to the effect of climate on growth
#' @param gmcsMinAge minimum age for which to predict full effect of growth/mortality -
#'                   younger ages are weighted toward a null effect with decreasing age
#' @param time the time of simulation - only used for writing outputs if applicable
#' @param cohortDefinitionCols cohortData columns that determine individual cohorts
#' @importFrom data.table `:=` as.data.table data.table fwrite rbindlist set setkey setnames 
#' @importFrom LandR asInteger
#' @importFrom terra ncell rast
#' @importFrom stats na.omit predict median
#' @rdname calculateClimateEffect
#' @export
calculateClimateEffect <- function(cohortData, pixelGroupMap, cceArgs,
                                   gmcsGrowthLimits, gmcsMinAge, time = NULL,
                                   cohortDefinitionCols = c("age", "speciesCode", "pixelGroup")) {
  
  originalCD <- cohortData
  cohortData <- copy(cohortData)
  # extract relevant args
  #assume that climate normals are present in cceArgs as cceArts$historicalClimate$historical_<var>_normal
  #assume climate variables are denoted by cceArgs$climateVariables
  if (!is.null(cceArgs$climateYear)) {
    thisYear <- cceArgs$climateYear
  } else {
    thisYear <- time
  }
  
  climateVariables <- cceArgs$climateVariablesForGMCS
  if (!is.null(cceArgs$currentClimateRasters)) {
    climateRasters <- cceArgs$currentClimateRasters[[climateVariables]]
  } else if (!is.null(cceArgs$projectedClimateRasters)) {
    climateRasters <- cceArgs$projectedClimateRasters
    climateRasters <- lapply(climateRasters[climateVariables], "[[", paste0("year", thisYear)) |>
      rast()
  } else {
    stop("no method of obtaining climate data for LandR.CS")
  }
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
    #TODO: think of ways this might get screwed up - prevent them 
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
  cohortData <- cohortData[age >= gmcsMinAge, .(speciesCode, B, mortality, age, pixelGroup)]
  
  #calculate biomass-weighted age, then log it
  cohortData[, standBiomass := sum(B), .(pixelGroup)]
  cohortData[, weightedAge := sum(B * age)/standBiomass, .(pixelGroup)]
  cohortData[, logAge := log(weightedAge)]
  
  #collapse cohorts by species
  cd <- cohortData[, .(biomass = sum(B)), .(pixelGroup, speciesCode, logAge, standBiomass)]
  
  #Make predictions
  predictedGrowth <- lapply(list(climateCovariates, historicalClimate), 
                            getModelPred, cohortData = cd, 
                            gmcsModel = cceArgs$gcsModel, 
                            stat = "growth")
  #take the 'historical' prediction for calculation of modifier
  subsetDT <- predictedGrowth[[2]][, .(pixelGroup, speciesCode, growth)]
  setnames(subsetDT, "growth", "pred_historicalGrowth")
  setnames(predictedGrowth[[1]], old = "growth", new = "pred_currentGrowth")
  #take means before join to avoid having to track fold
  
  predictedGrowth <- predictedGrowth[[1]][subsetDT, on = c("pixelGroup", "speciesCode")]
  predictedGrowth[, growthPred := asInteger(pred_currentGrowth/pred_historicalGrowth * 100)]
  
  #back to cohorts separated by age
  cohortData <- cohortData[predictedGrowth, on = c("pixelGroup", "speciesCode")]
  predictedMortality <- lapply(list(climateCovariates, historicalClimate), 
                               getModelPred, cohortData = cd, 
                               gmcsModel = cceArgs$mcsModel, 
                               stat = "mortality")
  #take the 'historical' prediction for calculation of modifier
  subsetDT <- predictedMortality[[2]][, .(pixelGroup, speciesCode, mortality)]
  
  #note cohort data has existing mortality column
  setnames(subsetDT, "mortality", "pred_historicalMortality")
  setnames(predictedMortality[[1]], old = "mortality", new = "pred_currentMortality")
  #take means before join to avoid having to track fold
  
  predictedMortality <- predictedMortality[[1]][subsetDT, on = c("pixelGroup", "speciesCode")]
  #mortality modifier is in g/m2
  #apply modifiers here 
  predictedMortality[pred_currentMortality < 0, pred_currentMortality := 0]
  predictedMortality[pred_historicalMortality < 0, pred_historicalMortality := 0]
  
  predictedMortality[, mortPred := asInteger(pred_currentMortality - pred_historicalMortality)]
  #need pooled biomass
  #this is the species-level biomass - must be called this for predict
  predictedMortality <- predictedMortality[cd[, .(pixelGroup, speciesCode, biomass)],
                                           on = c("pixelGroup", "speciesCode")]
  setnames(predictedMortality, old = "biomass", new = "speciesB")

  cohortData <- cohortData[predictedMortality, on = c("pixelGroup", "speciesCode")]
  
  cohortData[, growthPred := pmax(min(gmcsGrowthLimits), pmin(growthPred, max(gmcsGrowthLimits)))]
  #distribute mortality according to B 
  # B here is summed by species - so rename

  cohortData <- cohortData[, .(pixelGroup, speciesCode, age,  B, speciesB, standBiomass, growthPred, mortPred)]

  tooYoung <- originalCD[age < gmcsMinAge,]
  
  # speciesB does not include cohorts that are below the minimum age
  # this scales the prediction (in g/m2) by cohort biomass 
  # (for each species, for when 2 or more cohorts are present)
  cohortData[, mortPred := c(mortPred * B/speciesB)]
  # cohortData[, c("speciesB", "standBiomass") := NULL]
  
  if (is.null(cceArgs$climateYear)) {
    cohortData[, climateYear := time]
  } else {
    cohortData[, climateYear := cceArgs$climateYear]
  }
  #due to creation of summary module, decided to add these columns back in
  cohortData <- data.table::merge.data.table(cohortData, 
                    predictedGrowth[, .(pred_currentGrowth, pred_historicalGrowth, pixelGroup, speciesCode)])
  cohortData <- data.table::merge.data.table(cohortData, 
                                             predictedMortality[, .(pred_currentMortality, pred_historicalMortality, 
                                                                    pixelGroup, speciesCode)])
  cohortData <- climateCovariates[cohortData, on = c("pixelGroup")]
  #now add back in those that are too young
  #TODO first subset the columns of tooYoung so it doesnt' add bPM, bAP, mAge, etc. 
  tooYoung <- tooYoung[, intersect(names(tooYoung), names(cohortData)), with = FALSE]
  cohortData <- rbind(cohortData, tooYoung, fill = TRUE)
  cohortData[is.na(growthPred), c("growthPred", "mortPred") := .(100, 0)]
  
  return(cohortData = cohortData)
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
#' @importFrom data.table as.data.table copy melt.data.table rbindlist setnames
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