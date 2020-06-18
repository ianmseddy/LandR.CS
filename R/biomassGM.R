#'  CalculateClimateEffect
#'
#'  Predict biomass change with climate variables
#'
#' @param cceArgs a list of datasets used by the climate function
#' @param cohortData The LandR cohortData object
#' @param pixelGroupMap the pixelGroupMap needed to match cohorts with raster values
#' @param gmcsGrowthLimits lower and upper limits to the effect of climate on growth
#' @param gmcsMortLimits lower and upper limits to the effect of climate on mortality
#' @param gmcsMinAge mininmum age for which to predict growth/mortality
#' @importFrom data.table setkey data.table
#' @importFrom stats median
#' @importFrom raster getValues projection
#' @rdname calculateClimateEffect
#' @export
calculateClimateEffect <- function(cohortData, pixelGroupMap, cceArgs,
                                   gmcsGrowthLimits, gmcsMortLimits, gmcsMinAge){
  cohortData <- copy(cohortData)
  neededCols <- c("pixelGroup", 'speciesCode', 'age', 'B') %>%
    .[. %in% colnames(cohortData)]
  climCohortData <- cohortData[, ..neededCols]

  #extract relevant args
  ATA <- cceArgs$ATA
  CMI <- cceArgs$CMI
  CMInormal <- cceArgs$CMInormal
  mcsModel <- cceArgs$mcsModel
  gcsModel <- cceArgs$gcsModel

  if (ncell(CMI) != ncell(CMInormal)) {
    stop("different number of pixels in the climate data. Please review how these are created")
  }

  if (projection(ATA) != projection(CMInormal)) {
    stop("CRS of climate data is not identical. Please review how these are created")
  }

  CMIvals <- getValues(CMI)
  CMInormalvals <- getValues(CMInormal)
  ATAvals <- getValues(ATA)
  pixels <- getValues(pixelGroupMap)

  #Center observations on mean of original model data
  climateMatch <- data.table("pixelGroup" = pixels,
                             "CMI" = CMIvals,
                             "ATA" = ATAvals,
                             'CMInormal' = CMInormalvals)

  climateMatch <- climateMatch[!is.na(pixelGroup)] #Not all pixelGroups are in pixelGroupMap, because climCohortData is a subset
  #Take the median climate for each pixel group as some pixelgroups occur across multiple climate raster pixels
  climValues <- climateMatch[, .("CMI" = median(CMI, na.rm = TRUE),
                                 "ATA" = median(ATA, na.rm = TRUE),
                                 "CMInormal" = median(CMInormal, na.rm = TRUE)), by = "pixelGroup"]

  climCohortData[, logAge := log(age)]
  #set age = 0 to 1, to prevent -inf in prediction - this shouldn't affect predictions due to minimum age
  climCohortData[age == 0, logAge := 0]
  setkey(climCohortData, pixelGroup)
  setkey(climValues, pixelGroup)

  #Join cohort Data with climate data
  predData <- climCohortData[climValues]

  #remove NA values that exist only because of pixelGroupMap
  predData <- na.omit(predData)

  pixelGroupsPostSubset <- predData$pixelGroup
  agePostSubset <- predData$age
  speciesCodePostSubset <- predData$speciesCode

  predData <- predData[, .(logAge, ATA, CMI, CMInormal)]

  #Create the 'reference climate' dataset to normalize the prediction
  refClim <- predData
  refClim$CMI <- refClim$CMInormal #replace CMI with the CMI normal for 1950-2010
  refClim$ATA <- 0 #the anomaly by definition has 0 as nromal

  refClim[, CMInormal := NULL] #or the mortality model will be upset
  predData[, CMInormal := NULL]

  #make growth prediction as ratio
  growthPred <- asInteger(predict(gcsModel, predData, level = 0, asList = TRUE, type = "response")/
                            predict(gcsModel, refClim, level = 0, asList = TRUE, type = "response") * 100)
  growthPred[growthPred < min(gmcsGrowthLimits)] <- min(gmcsGrowthLimits)
  growthPred[growthPred > max(gmcsGrowthLimits)] <- max(gmcsGrowthLimits)

  #make mortality prediction
  mortPred <- asInteger(predict(object = mcsModel, parameter ='mu',
                                newdata = predData, level = 0, asList = TRUE, type = "response")/
                          predict(object = mcsModel, parameter = 'mu', newdata = refClim,
                                  level = 0, asList = TRUE, type = "response") * 100)

  mortPred[mortPred < min(gmcsMortLimits)] <- min(gmcsMortLimits)
  mortPred[mortPred > max(gmcsMortLimits)] <- max(gmcsMortLimits)


  if (anyNA(c(mortPred, growthPred))) {
    stop("error in climate prediction. NA value returned - this will break LANDR downstream")
  }

  climateEffect <- data.table("pixelGroup" = pixelGroupsPostSubset,
                              'speciesCode' = speciesCodePostSubset,
                              "age" = agePostSubset,
                              "growthPred" = growthPred,
                              "mortPred" = mortPred)

  if (length(cceArgs) > 5) {

    if (!any(is.null(cceArgs$transferTable), is.null(cceArgs$BECkey),
             is.null(cceArgs$currentBEC), is.null(cceArgs$ecoregionMap))) {
      #we do not want all the columns in cohortData, but we need ecoregionGroup and Provenance if present
      if (is.null(cohortData$Provenance)) {
        modCohortData <- cohortData[, .(pixelGroup, speciesCode, age, ecoregionGroup)]
      } else {
      modCohortData <- cohortData[, .(pixelGroup, speciesCode, age, ecoregionGroup, Provenance)]
      }

      geneticEffect <- calculateGeneticEffect(cohortData = modCohortData,
                                              BECkey = cceArgs$BECkey,
                                              pixelGroupMap = pixelGroupMap,
                                              transferTable = cceArgs$transferTable,
                                              ecoregionMap = cceArgs$ecoregionMap,
                                              currentBEC = cceArgs$currentBEC)

      setkey(climateEffect, pixelGroup, speciesCode, age)
      setkey(geneticEffect, pixelGroup, speciesCode, age)
      climateEffect <- geneticEffect[climateEffect]
      climateEffect[, growthPred := asInteger(HTp_pred * growthPred)]
      climateEffect[, HTp_pred := NULL] #get rid of this column

    } else {
      stop("cceArgs does not match methods available in LandR.CS")
    }
  }

  #restrict predictions to those above min stand age
  climateEffect[age < gmcsMinAge, mortPred := 100]
  climateEffect[age < gmcsMinAge, growthPred := 100]

  climateEffect <- climateEffect[cohortData[, .(pixelGroup, speciesCode, age)], on = c('pixelGroup', 'speciesCode', 'age')]
  #this is to fix any pixelGroups that were dropped by the na.omit of climData due to NA climate values

  climateEffect[is.na(growthPred), c('growthPred', 'mortPred') := .(100, 100)]

  #Because the params are numeric (e.g 66.667, the comparison forces the int to numeric)
  climateEffect[, c('growthPred', 'mortPred') := .(asInteger(growthPred), asInteger(mortPred))]

  return(climateEffect)
}


#'  own
#'  for predicting from gamlss with no random effect
#' @param fixed the fixed terms
#' @param random the random terms
#' @param correlation this is the correlation structure?
#' @param method Ceres help me
#' @param level the marginal or conditional predictor
#' @importFrom nlme lmeControl
#' @importFrom reproducible .grepSysCalls
#' @rdname own
#' @export
own <-function(fixed=~1, random = NULL, correlation = NULL, method = "ML",
               level = NULL, ...)
{
  #------------------------------------------
  # function starts here
  #------------------------------------------
  scall <- deparse(sys.call(), width.cutoff = 500L) #
  if (!is(fixed, "formula")) stop("fixed argument in lme() needs a formula starting with ~")
  #if (!is(random, "formula")) stop("formula argument in lme() needs a formula starting with ~")
  # we have to do somehing with corelation
  # if (!is.null(correlation)) {
  #   cor.for <- attr(correlation, "formula")
  #   if (!is.null(cor.for))
  #     cor.vars <- all.vars(cor.for)
  # }
  # else cor.vars <- NULL
  # get where "gamlss" is in system call
  # it can be in gamlss() or predict.gamlss()

  # use `tail` because there may be more than one gamlss, e.g., with system.call(gamlss).
  #  take the last one ... thus tail(..., 1)
  position <- tail(reproducible::.grepSysCalls(sys.calls(), "gamlss"), 1)
  # for (i in length(rexpr):1)
  # {
  #   position <- i # get the position
  #   if (rexpr[i]==TRUE) break
  # }
  #
  gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
  ##--- get the lme control values
  control  <- lmeControl(...)
  ## get the data
  if (sys.call(position)[1]=="predict.gamlss()")
  { # if predict is used
    Data <- get("data", envir=gamlss.env)
  }
  else if (sys.call(position)[1]=="gamlss()")
  { # if gamlss() is used
    if (is.null(get("gamlsscall", envir=gamlss.env)$data))
    { # if no data argument but the formula can be interpreted
      Data <- model.frame(formula)
    }
    else
    {# data argument in gamlss
      Data <- get("gamlsscall", envir=gamlss.env)$data
    }
  }
  else  {
    Data <- get("data", envir=gamlss.env)
  }
  Data <- if (any(attributes(eval(substitute(Data)))$class=="groupedData")) eval(substitute(Data))
  else data.frame(eval(substitute(Data)))
  #=====
  len <- dim(Data)[1] # get the lenth of the data
  ## out
  xvar <- rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
  attr(xvar,"fixed")       <- fixed
  attr(xvar,"random")      <- random
  attr(xvar,"method")      <- method
  attr(xvar,"correlation") <- correlation
  attr(xvar,"level")       <- level
  attr(xvar,"control")     <- control
  attr(xvar, "gamlss.env") <- gamlss.env
  if (any(attributes(Data)$class=="groupedData")) {
    attr(xvar, "data") <- Data } else {
      attr(xvar, "data") <- as.data.frame(Data)
    }
  attr(xvar, "call")       <- substitute(gamlss.own(data[[scall]], z, w, ...))
  attr(xvar, "class")      <- "smooth"
  xvar
}



#' gamlss.own
#' the definition of the backfitting additive function, Authors: Mikis Stasinopoulos, Marco Enea
#' @param x description missing
#' @param y description missing
#' @param w description missing
#' @importFrom nlme lme
#' @rdname gamlss.own
#' @export
gamlss.own <- function(x, y, w, xeval = NULL, ...)
{
  fixed <- attr(x, "fixed")
  random <- attr(x, "random")
  correlation <- attr(x, "correlation")
  method <- attr(x, "method")
  level <- attr(x, "level")
  fix.formula <-
    as.formula(paste("Y.var", deparse(fixed, width.cutoff = 500L), sep = ""))
  control <- as.list(attr(x, "control"))
  #gamlss.env <- as.environment(attr(x, "gamlss.env"))
  OData <- attr(x, "data")
  Data <-  if (is.null(xeval))
    OData #the trick is for prediction
  else
    OData[seq(1, length(y)), ]
  if (any(attributes(Data)$class == "groupedData")) {
    Data$W.var <- 1 / w
    Data$Y.var <- y
  } else {
    Y.var <- y
    W.var <- 1 / w
    Data <- data.frame(eval(substitute(Data)), Y.var, W.var)
  }
  #       Data <- data.frame(eval(substitute(Data)),y,wei=1/w)
  # fit  <-  lme(All$fixed, data = Data, random=All$random, weights=varFixed(~wei),  method="ML")
  #
  #              (fixed, data = sys.frame(sys.parent()), random, correlation = NULL,
  #           weights = NULL, subset, method = c("REML", "ML"), na.action = na.fail,
  #           control = list(), contrasts = NULL, keep.data = TRUE)
  # lme(fixed = fixed, random = random, data = data,
  #    correlation = correlation, control = control, weights = varFixed(w.formula),
  #    method = "ML", ...)
  fit <- lme(fixed = fix.formula,
             data = Data,
             random = random,
             weights = varFixed( ~ W.var),
             correlation = correlation,
             control = control,
             method = method
  )
  fv <- fitted(fit)
  residuals <- y - fv
  N <- sum(w != 0)
  df <-  N - (sum(w * (y - fv) ^ 2)) / (fit$sigma ^ 2)
  if (is.null(xeval)){
    list(
      fitted.values = fv,
      residuals = residuals,
      nl.df = df - 1,
      lambda = fit$sigma,
      # Mikis 10-6-19 df should be df-1 not df
      coefSmo = fit,
      var = NA
    )    # var=fv has to fixed
  }
  else {
    # ll<-dim(OData)[1]
    # assign("fix.formula",fix.formula,envir=globalenv())
    # on.exit(rm(fix.formula,envir=globalenv()))
    # pred <- eval(expression(predict(fit,newdata = OData[seq(length(y)+1,ll),])),envir=environment() )
    fit$call$fixed <- substitute(fix.formula)
    ll <- dim(OData)[1]
    pred <-
      if (is.null(level))
        predict(fit, newdata = OData[seq(length(y) + 1, ll), ])
    else
      predict(fit, newdata = OData[seq(length(y) + 1, ll), ], level = level)
  }
}


#'  calculateGeneticEffect
#'
#'  Predict climate induced height reduction due to genetics
#'
#' @param BECkey a key that matches BECraster code with transfer table
#' @param cohortData The LandR cohortData object
#' @param pixelGroupMap the pixelGroupMap needed to match cohorts with raster values
#' @param transferTable a table with genetic performance of species in each variant
#' @param currentBEC the current projected BEC zone
#' @param ecoregionMap a raster an RAT that matches ecoregionGroup to ecoregion
#' @importFrom data.table setkey data.table
#' @importFrom stats median
#' @importFrom raster getValues projection
#' @rdname calculateGeneticEffect
#' @export
calculateGeneticEffect <- function(BECkey, cohortData, pixelGroupMap, transferTable, currentBEC, ecoregionMap){

  transferTable <- copy(transferTable) #this is necessary due to column name changes
  BECkey <- copy(BECkey) #this is necessary due to class change
  #1. get BEC zones of each ecoregionGroup
  ecoregionKey <- as.data.table(ecoregionMap@data@attributes[[1]])
  setnames(ecoregionKey, 'ID', 'ecoregionMapCode') #Change ID, because ID in BECkey = ecoregion, not mapcode
  BECkey[, ID := as.factor(as.character(ID))]

  ecoregionKey <- BECkey[ecoregionKey, on = c("ID" = 'ecoregion')] #now we have zsv of cohortData$ecoregionGroup
  ecoregionKeySmall <- ecoregionKey[, .(zsv, ecoregionGroup)]

#2. Find Provenance of cohortData
  bugCatch <- nrow(cohortData)
  if (is.null(cohortData$Provenance)){
    cohortData <- cohortData[, .(speciesCode, ecoregionGroup, pixelGroup, age)] %>%
      ecoregionKeySmall[., on = c("ecoregionGroup" = 'ecoregionGroup')]
    setnames(cohortData, 'zsv', 'Provenance')
  } else {
    cohortData <- cohortData[, .(speciesCode, ecoregionGroup, pixelGroup, age, Provenance)] %>%
      ecoregionKeySmall[., on = c("ecoregionGroup" = 'ecoregionGroup')]
    setnames(cohortData, 'zsv', 'assumedProvenance')
    cohortData[is.na(Provenance), Provenance := assumedProvenance]
    cohortData[, assumedProvenance := NULL] #Confirm this doesn't erase parts of Provenance by reference
  }

#3. Assign the mode among projected BECs for each pixelGroup
  projBEC <- data.table(pixelGroup = getValues(pixelGroupMap), BEC = getValues(currentBEC)) %>%
    na.omit(.) %>%
    .[, BEC := as.factor(BEC)]
  projBEC <- projBEC[, .N, .(pixelGroup, BEC)]
  projBEC[, modeBEC := max(N), .(pixelGroup)]
  projBEC <- projBEC[N == modeBEC] %>%
    .[, c('modeBEC', 'N') := NULL]

  #Find ties
  counts <- projBEC[, .N, .(pixelGroup)]
  noTies <- projBEC[pixelGroup %in% counts[N == 1]$pixelGroup]
  ties <- projBEC[pixelGroup %in% counts[N > 1]$pixelGroup]
  rm(counts)
  #randomly order, then remove duplicates
  if (nrow(ties) > 1) {
  ties$foo <- sample(x = 1:nrow(ties), size = nrow(ties))
  setkey(ties, foo)
  ties <- ties[!duplicated(ties[, .(pixelGroup)])]
  ties[, foo := NULL]
  assignedBEC <- rbind(ties, noTies)
  } else {
    assignedBEC <- noTies
  }

  if (nrow(assignedBEC) != length(unique(projBEC$pixelGroup))) {
    stop("Error: mismatch in pixelGroups and projected BECs, debug LandR.CS")
  }
  rm(ties, noTies, projBEC)

#4.Join tables
  assignedBEC <- BECkey[assignedBEC, on = c("ID" = 'BEC')] %>%
    .[, .(pixelGroup, zsv)]
  setnames(assignedBEC, old = 'zsv', new = 'currentClimate')
  cohortData <- assignedBEC[cohortData, on = c('pixelGroup' = 'pixelGroup')]

#5.Add Provenance if missing
  #I think this is redundant
  cohortData <- ecoregionKeySmall[cohortData, on = c("ecoregionGroup" = 'ecoregionGroup')]
  if (is.null(cohortData$Provenance)) {
    cohortData[, Provenance := zsv]
  } else {
    cohortData[is.na(Provenance), Provenance := zsv]
  }
  cohortData[, zsv := NULL]

  setnames(transferTable, old = c("BECvarfut_plantation", 'BECvar_seed'), new = c("currentClimate", "Provenance"))

  cohortData <- transferTable[cohortData, on = c('currentClimate' = 'currentClimate',
                                                  'Provenance' = 'Provenance',
                                                  'speciesCode' = 'speciesCode')] %>%
    .[, .(pixelGroup, speciesCode, age, Provenance, HTp_pred)]

  if (any(is.na(cohortData$Provenance))) {
    browser()
    #what the hell happened?
  }

  if (nrow(cohortData) != bugCatch) {
    stop("unequal row count after calculate genetic effect. debug LandR.CS")
  }
  return(cohortData)
}
