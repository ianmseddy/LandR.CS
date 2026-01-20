globalVariables(c(
  ".", "age", "logAge", "pixelGroup", "..cohortDefinitionCols", "..addedColumns", "..neededCols", ":="
  "spp", "speciesCode", "obs", "pred", "resid"
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

  #protect yourself at all times
  cohortData <- copy(cohortData)
   # extract relevant args
  browser()
  #assume that climate normals are present in cceArgs as cceArts$historicalClimate$historical_<var>_normal
  #assume climate variables are denoted by cceArgs$climateVariables
  climateVariables <- cceArgs$climateVariablesForGMCS
  climateRasters <- cceArgs$projectedClimateRasters
   rasterToGet <- paste0("year", cceArgs$climateYear)
  #make composite raster
  climateRasters <- lapply(climateRasters[climateVariables], "[[", rasterToGet) |>
    rast()
  
  #get the variables which require anomalies
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
  #there is probably a clever way 
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
                    .SDcol = names(allClimateRas), .(pixelGroup)]
  
  #don't forget to.....
  #TODO: this is WRONG - calculate the biomass-weighted stand age, take the log of that. Identical for a PG
  #TODO: also sum biomass by species within a PG 
  
  #subset cd for necessary columns
  cohortData <- cohortData[, .(speciesCode, B, mortality, age, pixelGroup)]
  
  #calculate biomass-weighted age- then log it
  cohortData[, standBiomass := sum(B), .(pixelGroup)]
  #TODO: confirm this is correct when you have uneven-aged stands
  cohortData[, weightedAge := sum(B * age)/sumB, .(pixelGroup)]
  cohortData[, logAge := log(weightedAge)]

  #we can predict species but not individual cohorts (this is absent from PSP, as we do not have their ages)
  cohortData <- cohortData[, .(biomass = sum(B)), .(pixelGroup, speciesCode, logAge, standBiomass)]
  
  predData <- cohortData[climateCovariates, on = "pixelGroup"]
  setnames(predData, "speciesCode", "spp") #to match fitting
  #Next up - take the character vector and make into dummy variables like so
  #give each species "spp" at the start 
  #you may need to create the extra dummy columns for species that aren't present, e.g. sppTsug_het
  predData <- model.matrix(~ . + 0, data = predData)
  predData <- as.data.table(predData)
  
  #if tree species are absent, they must be added for the predict
  fittedCols <- names(cceArgs$gcsModel$Fold1$valData)
  missing <- setdiff(fittedCols, c(names(predData), "obs", "pred", "resid"))
  if (length(missing)){
    set(predData, j = missing, value =  0)
  }
  
  growthPreds <- lapply(cceArgs$gcsModel, function(Fold, data = predData){
    out <- copy(data)
    out[, predictedGrowth := predict(object = Fold$mod, newdata = data)]
    return(out)
  }) |>
    rbindlist()
  #now summarize 
  #probably you will have to do this 4 times so turn it into a function
  out <- melt.data.table(growthPreds, id.vars = c("pixelGroup", "predictedGrowth"), 
                  measure.vars = grep("spp", colnames(growthPreds), value = TRUE),
                  variable.name = "speciesCode", value.name = "presence")
  out <- out[presence == 1,] #keep only presences
  out[, speciesCode := gsub(pattern = "spp", replacement = "", speciesCode)]
  out[, speciesCode := factor(speciesCode, levels = levels(cohortData$speciesCode))]
  
  
  ####OLD###
  # if (length(anomalies) > 0) {
  #   #check if year is already subset - if not, subset to PSP period
  #   anomalyVar <-
  #   setnames(anomalyData, old = anomalies, new = names(anomalies))
  #   PSPmodelData <- anomalyData[PSPmodelData, on = c("OrigPlotID1")]
  # 
  #   #recalculate the anomaly(s) via subtraction
  #   for (i in 1:length(anomalies)) {
  #     #index because the name isn't preserved if you take the object itself
  #     temp <- anomalies[i]
  #     setnames(PSPmodelData, c(temp, names(temp)), c("var", "anom"))
  #     PSPmodelData[, anom := var - anom]
  #     newOrderCols <- names(PSPmodelData)[!names(PSPmodelData) %in% c("var", "anom") ]
  #     setcolorder(PSPmodelData, newOrderCols)
  #     setnames(PSPmodelData, c("var", "anom"), c(temp, names(temp)))
  #   }
  # }
  # 
  # 
  # 
  # 
  # 
  # ATA <- cceArgs$projectedClimateRasters[["ATA"]][[paste0("year", year)]]
  # CMI <- cceArgs$projectedClimateRasters[["CMI"]][[paste0("year", year)]]
  # CMInormal <- cceArgs$historicalClimateRasters[["CMI_normal"]]
  # mcsModel <- cceArgs$mcsModel
  # gcsModel <- cceArgs$gcsModel
  # 
  # cohortData <- Copy(cohortData)
  # neededCols <- c(cohortDefinitionCols, "B")
  # neededCols <- neededCols[neededCols %in% colnames(cohortData)]
  # climCohortData <- cohortData[, ..neededCols]
  # 
  # if (!compareGeom(CMI, CMInormal, ATA, stopOnError = FALSE)) {
  #   stop("different number of pixels in the climate data. Please review how these are created")
  # }
  # 
  # 
  # CMIvals <- as.vector(CMI)
  # CMInormalvals <- as.vector(CMInormal)
  # ATAvals <- as.vector(ATA)
  # pixels <- as.vector(pixelGroupMap)
  # 
  # # Center observations on mean of original model data
  # climateMatch <- data.table(
  #   "pixelGroup" = pixels,
  #   "CMI" = CMIvals,
  #   "ATA" = ATAvals,
  #   "CMInormal" = CMInormalvals
  # )
  # 
  # climateMatch <- climateMatch[!is.na(pixelGroup)]
  # # Not all pixelGroups are in pixelGroupMap, because climCohortData is a subset
  # # Take the median climate for each pixel group as some pixelgroups occur across multiple climate raster pixels
  # climValues <- climateMatch[, .(
  #   "CMI" = median(CMI, na.rm = TRUE),
  #   "ATA" = median(ATA, na.rm = TRUE),
  #   "CMInormal" = median(CMInormal, na.rm = TRUE)
  # ), by = "pixelGroup"]
  # 
  # climCohortData[, logAge := log(age)]
  # # set age = 0 to 1, to prevent -inf in prediction - this shouldn't affect predictions due to minimum age
  # climCohortData[age == 0, logAge := 0]
  # setkey(climCohortData, pixelGroup)
  # setkey(climValues, pixelGroup)
  # 
  # # Join cohort Data with climate data
  # predData <- climCohortData[climValues]
  # 
  # # remove NA values that exist only because of pixelGroupMap
  # predData <- na.omit(predData)
  # 
  # pixelGroupsPostSubset <- predData$pixelGroup
  # agePostSubset <- predData$age
  # speciesCodePostSubset <- predData$speciesCode
  # 
  # modCohortDef <- FALSE
  # 
  # # necessary for joining if cohortData has added columns
  # if (length(cohortDefinitionCols[!cohortDefinitionCols %in% c("age", "pixelGroup", "speciesCode")]) > 0) {
  #   modCohortDef <- TRUE
  #   addedColumns <- cohortDefinitionCols[!cohortDefinitionCols %in% c("age", "pixelGroup", "speciesCode")]
  #   addedColumns <- predData[, ..addedColumns]
  # }
  # 
  # predData <- predData[, .(logAge, ATA, CMI, CMInormal)]
  # 
  # # Create the 'reference climate' dataset to normalize the prediction
  # refClim <- predData
  # refClim$CMI <- refClim$CMInormal # replace CMI with the CMI normal for 1950-2010
  # refClim$ATA <- 0 # the anomaly by definition has 0 as nromal
  # 
  # refClim[, CMInormal := NULL] # or the mortality model will be upset
  # predData[, CMInormal := NULL]
  # 
  # # make growth prediction as ratio
  # growthPred <- asInteger(predict(gcsModel, predData, level = 0, asList = TRUE, type = "response") /
  #                           predict(gcsModel, refClim, level = 0, asList = TRUE, type = "response") * 100)
  # growthPred[growthPred < min(gmcsGrowthLimits)] <- min(gmcsGrowthLimits)
  # growthPred[growthPred > max(gmcsGrowthLimits)] <- max(gmcsGrowthLimits)
  # growthPred[is.na(growthPred)] <- max(gmcsGrowthLimits)
  # 
  # # make mortality prediction
  # 
  # mortPred <- asInteger(predict(
  #   object = mcsModel, parameter = "mu",
  #   newdata = predData, level = 0, asList = TRUE, type = "response"
  # ) /
  #   predict(
  #     object = mcsModel, parameter = "mu", newdata = refClim,
  #     level = 0, asList = TRUE, type = "response"
  #   ) * 100)
  # mortPred[mortPred < min(gmcsMortLimits)] <- min(gmcsMortLimits)
  # mortPred[mortPred > max(gmcsMortLimits)] <- max(gmcsMortLimits)
  # 
  # if (anyNA(c(mortPred, growthPred))) {
  #   mortPred[is.na(mortPred)] <- max(gmcsMortLimits)
  # 
  #   warning("NA in climate prediction. Likely integer overflow - setting to gmcsLimits")
  # }
  # 
  # # predict requires exact same columns in data.frame at the moment, hence this clumsy rebuilding
  # climateEffect <- data.table(
  #   "pixelGroup" = pixelGroupsPostSubset,
  #   "speciesCode" = speciesCodePostSubset,
  #   "age" = agePostSubset,
  #   "growthPred" = growthPred,
  #   "mortPred" = mortPred
  # )
  # if (modCohortDef) {
  #   climateEffect <- cbind(climateEffect, addedColumns)
  # }
  # 
  # # restrict predictions to those above min stand age
  # climateEffect[age < gmcsMinAge, growthPred := as.integer(100 + ((growthPred - 100) * (age / gmcsMinAge)))]
  # climateEffect[age < gmcsMinAge, mortPred := as.integer(100 + ((mortPred - 100) * (age / gmcsMinAge)))]
  # temp <- cohortData[, ..cohortDefinitionCols]
  # climateEffect <- climateEffect[temp, on = cohortDefinitionCols]
  # rm(temp)
  # # this is to fix any pixelGroups that were dropped by the na.omit of climData due to NA climate values
  # # which should be quite rare but persist with postProcess problems
  # climateEffect[is.na(growthPred), c("growthPred", "mortPred") := .(100, 100)]
  # 
  # # Because the params are numeric (e.g 66.667, the comparison forces the int to numeric)
  # climateEffect[, c("growthPred", "mortPred") := .(asInteger(growthPred), asInteger(mortPred))]

  return(climateEffect)
}

#' getModelPred
#'
#' Join cohortData with either the growth or mortality xgboost models
#' This is a utility function to avoid repetition as 2 predictions are needed
#' for each model (under the current climate, and historical) to achieve the 
#' climate modifier

#' @param cohortData The LandR cohortData object
#' @param climateCovariates data.table of mean climate covariates by pixelGroup 
#' @param cohortDefinitionCols cohortData columns that determine individual cohorts
#' @importFrom data.table setkey data.table
#' @importFrom reproducible Copy
#' @importFrom LandR asInteger
#' @importFrom terra compareGeom ncell
#' @importFrom data.table melt.data.table copy
#' @rdname getModelPred
getModelPred <- function(cohortData, climateCovariates, gmcsModel) {
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
}
