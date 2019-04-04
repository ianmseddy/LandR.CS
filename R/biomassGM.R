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
  if (is.null(CMD)) {
    stop("Missing climate data needed to run LandR.CS - consider running module gmcsDataPrep and PSP_Clean")
  }
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
  growthPred <- predict(gcsModel, predData,
                        level = 0, asList = TRUE)
  if (anyNA(growthPred)) {
    browser()
  }
  #Prediction is in tons/ha, must be rescaled to g/m2
  #1000000 g/10000 m2 = 100 * g/m2
  growthPred <- growthPred * 100

  mortPred <- predict(mcsModel, predData, level = 0, asList = TRUE)
  if (anyNA(growthPred)) {
    browser()
  }
  mortPred <- mortPred * 100

  climateEffect <- data.table("pixelGroup" = predData$pixelGroup,
                              "growthPred" = growthPred,
                              "mortPred" = mortPred)
  #good work Ian
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
