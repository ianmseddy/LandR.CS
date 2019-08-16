#'  CalculateClimateEffect
#'
#'  Predict biomass change with climate variables
#'
#' @param cohortData The LandR cohortData object
#' @param CMI A raster with annual Cumulative Moisture Index values
#' @param ATA A raster with annual temperature anomaly values
#' @param gcsModel climate-sensitive growth mixed effect model object created by gmcsDataPrep
#' @param mcsModel climate-sensitive mortality mixed effect model object created by gmcsDataPrep
#' @param pixelGroupMap the pixelGroupMap needed to match cohorts with raster values
#' @param centeringVec the means of the data used to create gcsModel and mcsModel
#' @param CMInormal raster of CMI normals for 1950-2010
#' @importFrom data.table setkey data.table
#' @importFrom stats median
#' @importFrom raster getValues
#' @rdname calculateClimateEffect
#' @export
calculateClimateEffect <- function(cohortData, CMI, ATA, gcsModel, mcsModel,
                                   pixelGroupMap, centeringVec, CMInormal){
  if (is.null(CMI)) {
    stop("Missing climate data needed to run LandR.CS - consider running module gmcsDataPrep and PSP_Clean")
  }
  browser()

  CMIvals <- getValues(CMI)
  CMInormalvals <- getValues(CMInormal)
  ATAvals <- getValues(ATA)
  pixels <- getValues(pixelGroupMap)

  #if all values are 0, its because the current time is before 2011
  if (all(unique(CMIvals) == 0, na.rm = TRUE)) {
    return(NULL)
  }

  #Center observations on mean of original model data
  climateMatch <- data.table("pixelGroup" = pixels,
                             "CMI" = CMIvals,
                             "ATA" = ATAvals,
                             'CMInormal' = CMInormalvals)

  climateMatch <- climateMatch[!is.na(pixelGroup)]
  #Take the median climate variables for each pixel group
  out <- climateMatch[, list("CMI" = median(CMI, na.rm = TRUE),
                             "ATA" = median(ATA, na.rm = TRUE),
                             "CMInormal" = median(CMInormal)), by = "pixelGroup"]

  #summarize cohortData by biomass
  cohortData <- cohortData[, list(age = max(age), B = sum(B)), by = "pixelGroup"]
  cohortData$logAge <- log(cohortData$age)
  setkey(cohortData, pixelGroup)
  setkey(out, pixelGroup)
  #Join cohort Data with climate data
  predData <- out[cohortData]

  #Create the 'avg climate' dataset to normalize the prediction
  avgClim <- predData
  avgClim$CMI <- predData$CMInormal #replace CMI with the CMI normal for 1950-2010
  avgClim$ATA <- 0 #the anomaly by definition has 0 as nromal

  #make growth prediction
  growthPred <- predict(gcsModel, predData, level = 0, asList = TRUE) -
    predict(gcsModel, avgClim, level = 0, asList = TRUE)

  #make mortality prediction
  mortPred <- predict(mcsModel, predData, level = 0, asList = TRUE)
  #back transform
  mortPred <- exp(mortPred) - centeringVec['minMort']

  #predict the 'average climate' mortality
  nullMort <- predict(mcsModel, avgClim, level = 0, asList = TRUE)
  #back transform
  nullMort <- exp(nullMort) - centeringVec['minMort']

  mortPred <- mortPred - nullMort

  if (anyNA(mortPred)) {
    stop("error in climate mortality prediction. NA value returned")
  }

  climateEffect <- data.table("pixelGroup" = predData$pixelGroup,
                              "growthPred" = growthPred,
                              "mortPred" = mortPred)
  return(climateEffect)
}


#'  assignClimateEffect
#'
#'  Calculates climate-dependent mortality and growth proportional to
#'  each cohorts biomass in the pixelGroup
#'
#' @param subCohortData The LandR cohortData object
#' @param predObj climate prediction object
#' @param type predict growth or mortality
#' @export
#' @rdname assignClimateEffect
assignClimateEffect <- function(subCohortData, predObj, type){

  #if predObj is null, the time must be before 2011.
  if (is.null(predObj)) {
    return(0)
  }

  subCohorts <- subCohortData[,.("pixelGroup" = pixelGroup, "B" = B, 'aNPPAct' = aNPPAct)]

  #Mortality is proportional to each cohort's biomass
  if (type == 'mortPred') {
    subCohorts[, "sumB" := sum(B), by = "pixelGroup"]
    subCohorts[, "propB" := B/sumB,]
    #subCohortData should be sorted on pixelGroup. Need to preserve original order
    subCohorts[, "rowOrder" := as.numeric(row.names(subCohorts))]
    setkey(subCohorts, "pixelGroup")
    setkey(predObj, "pixelGroup")
    subCohorts <- predObj[subCohorts]
    subCohorts[, "climStat" := eval(parse(text = type)) * propB]
  } else {
    #Growth is proportionalt o each cohort's annual net primary productivity
    subCohorts[, "sumANPP" := sum(aNPPAct), by = "pixelGroup"]
    subCohorts[, "propANPP" := aNPPAct/sumANPP,]
    #subCohortData should be sorted on pixelGroup. Need to preserve original order
    subCohorts[, "rowOrder" := as.numeric(row.names(subCohorts))]
    setkey(subCohorts, "pixelGroup")
    setkey(predObj, "pixelGroup")
    subCohorts <- predObj[subCohorts]
    subCohorts[, "climStat" := eval(parse(text = type)) * propANPP]
  }
  setkey(subCohorts, "rowOrder") #Back to original order
  return(subCohorts$climStat)

}
