utils::globalVariables(c(
  ".SD", ".", "age", "aNPPAct", "B", "ClimateYear", "growth", "growthPred", 
  "pred_currentGrowth", "pred_currentMortality", "pred_historicalGrowth", "pred_historicalMortality", 
  "logAge",  "mortality", "mortPred", "obs", "pixelGroup", "pred", "predCol", "presence", 
  "spp", "speciesB", "speciesCode", "resid", "standBiomass", "weightedAge", "year"
))
# 
# .onLoad <- function(libname, pkgname) {
#   op <- options()
#   op.LandR.CS <- list(
#     LandR.CS.debug = FALSE,
#     LandR.CS.logPath = "outputs/LandRCS"
#   )
#   
#   toset <- !(names(op.LandR.CS) %in% names(op))
#   if (any(toset)) options(op.LandR.CS[toset])
# }