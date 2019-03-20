#'  CalculateClimateEffect
#'
#'  Predict biomass change with climate variables
#'
#' @param cohortData The LandR cohortData object
#' @param CMD A raster with annual Cumulative Moisture Deficit values
#' @param ATA A raster with annual temperature anomaly values
#' @param gcsModel climate-sensitive growth mixed effect model object created by gmcsDataPrep
#' @param mcsModel climate-sensitive mortality mixed effect model object created by gmcsDataPrep
#' @param pixelGroupMap the pixelGroupMap needed to match cohorts with raster values
#' @param centeringVec the means of the data used to create gcsModel and mcsModel
#' @importFrom data.table setkey data.table
#' @importFrom stats median
#' @importFrom raster getValues
#' @rdname calculateClimateEffect
#' @export
calculateClimateEffect <- function(cohortData, CMD, ATA, gcsModel, mcsModel,
                                   pixelGroupMap, centeringVec){

  CMDvals <- getValues(CMD)
  ATAvals <- getValues(ATA)
  pixels <- getValues(pixelGroupMap)
  #Center observations on mean of original model data
  climateMatch <- data.table("pixelGroup" = pixels,
                             "mCMD" = CMDvals - centeringVec["CMD"],
                             "mATA" = ATAvals - centeringVec["ATA"])

  climateMatch <- climateMatch[!is.na(pixelGroup)]
  #Take the median climate variables for each pixel group
  out <- climateMatch[, list("mCMD" = median(mCMD, na.rm = TRUE), "mATA" = median(mATA, na.rm = TRUE)), by = "pixelGroup"]

  #summarize cohortData by biomass
  cohortData <- cohortData[, list(age = max(age), B = sum(B)), by = "pixelGroup"]
  cohortData$mLogAge <- log(cohortData$age) - centeringVec["logAge"]
  setkey(cohortData, pixelGroup)
  setkey(out, pixelGroup)
  #Join cohort Data with climate data
  predData <- out[cohortData]

  #make prediction
  growthPred <- predict(gcsModel, predData, na.rm = TRUE,
                        level = 0, asList = TRUE)
  #Prediction is in tons/ha, must be rescaled to g/m2
  #1000000 g/10000 m2 = 100 * g/m2
  growthPred <- growthPred * 100

  mortPred <- predict(mcsModel, predData, na.rm = TRUE,
                      level = 0, asList = TRUE)
  mortPred <- mortPred * 100

  climateEffect <- data.table("pixelGroup" = predData$pixelGroup,
                              "growthPred" = growthPred,
                              "mortPred" = mortPred)
  #good work Ian
   return(climateEffect)
}

#'  Calculate climate growth
#' @param subCohortData The LandR cohortData object
#' @param predObj climate prediction object
#' @rdname calculateClimateMortality
#'
#' @export
calculateClimateMortality <- function(subCohortData, predObj){

  subCohorts <- subCohortData[,.("pixelGroup" = pixelGroup, "B" = B)]
  subCohorts[, "sumB" := sum(B), by = "pixelGroup"]
  subCohorts[, "propB" := B/sumB,]
  #subCohortData should be sorted on pixelGroup. Need to preserve original order
  subCohorts[, "rowOrder" := as.numeric(row.names(subCohorts))]
  setkey(subCohorts, "pixelGroup")
  setkey(predObj, "pixelGroup")
  predObj[subCohorts]
  subCohorts[, "climMortality" := mortPred * propB]
  setkey(subCohorts, "rowOrder") #Back to original order
  return(subCohorts$climMortality)
}

#'  CalculateClimateGrowth
#' @param subCohortData The LandR cohortData object
#' @param predObj climate prediction object
#' @export
#' @rdname calculateClimateGrowth
calculateClimateGrowth <- function(subCohortData, predObj){

  subCohorts <- subCohortData[,.("pixelGroup" = pixelGroup, "B" = B)]
  subCohorts[, "sumB" := sum(B), by = "pixelGroup"]
  subCohorts[, "propB" := B/sumB,]
  #subCohortData should be sorted on pixelGroup. Need to preserve original order
  subCohorts[, "rowOrder" := as.numeric(row.names(subCohorts))]
  setkey(subCohorts, "pixelGroup")
  setkey(predObj, "pixelGroup")
  predObj[subCohorts]
  subCohorts[, "climGrowth" := growthPred * propB]
  setkey(subCohorts, "rowOrder") #Back to original order
  return(subCohorts$climGrowth)

}
