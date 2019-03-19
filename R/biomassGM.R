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
#' @importFrom data.table setkey data.table
#' @importFrom stats median
#' @importFrom raster getValues
#' @rdname calculateClimateEffect
#' @export
calculateClimateEffect <- function(cohortData, CMD, ATA, gcsModel, mcsModel, pixelGroupMap){
  browser()
  CMDvals <- getValues(CMD)
  ATAvals <- getValues(ATA)
  pixels <- getValues(pixelGroupMap)
  #Center observations on mean
  climateMatch <- data.table("pixelGroup" = pixels,
                             "mCMD" = CMDvals - mean(CMDvals, na.rm = TRUE),
                             "mATA" = ATAvals - mean(ATAvals, na.rm = TRUE))

  climateMatch <- climateMatch[!is.na(pixelGroup)]
  out <- climateMatch[, list("mCMD" = median(mCMD, na.rm = TRUE), "mATA" = median(mATA, na.rm = TRUE)), by = "pixelGroup"]

  #summarize cohortData by biomass
  cohortData <- cohortData[, list(age = max(age), B = sum(B)), by = "pixelGroup"]
  cohortData$mLogAge <- log(cohortData$age) - mean(log(cohortData$age))
  setkey(cohortData, pixelGroup)
  setkey(out, pixelGroup)
  predData <- out[cohortData]
  #THE PREDICTION IS IN TONS PER HECTARE NEVER FORGET
  newPred <- predict(gmcsModel, predData,
                        na.rm = TRUE, level = 0, asList = TRUE)

  #Hurray now we have the prediction whoo hoooo. But is it still in units of t/ha, so multiply by 100.

   return(climateMatch)

  #Need to generate predicted changes in biomass, return that object.
  #Other 2 functions will return proportional changes in mortality and growth

}

#'  Calculate climate growth
#' @param cohortData The LandR cohortData object
#' @param predObj climate prediction object
#' @rdname calculateClimateMortality
#'
#' @export
calculateClimateMortality <- function(cohortData, predObj){
  print(predObj)
  cohortData #what needs to be returned is a vector
  #Join cohort data? Take max age, convert to log, take mean.

  #Need to generate predicted changes in biomass, return that object.
  #Other 2 functions will return proportional changes in mortality and growth

}

#'  CalculateClimateGrowth
#' @param cohortData The LandR cohortData object
#' @param predObj climate prediction object
#' @export
#' @rdname calculateClimateGrowth
calculateClimateGrowth <- function(cohortData, predObj){
  cohortData #what needs to be returned is a vector
  #Join cohort data? Take max age, convert to log, take mean.

  #Need to generate predicted changes in biomass, return that object.
  #Other 2 functions will return proportional changes in mortality and growth

}
