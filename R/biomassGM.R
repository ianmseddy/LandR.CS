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
#' @param CMInormal raster of CMI normals the reference period
#' @param gmcsGrowthLimits lower and upper limits to the effect of climate on growth
#' @param gmcsMortLimits lower and upper limits to the effect of climate on mortality
#' @param gmcsMinAge mininmum age for which to predict growth/mortality
#' @importFrom data.table setkey data.table
#' @importFrom stats median
#' @importFrom raster getValues
#' @importFrom LandR projection
#' @rdname calculateClimateEffect
#' @export
calculateClimateEffect <- function(cohortData, CMI, ATA, gcsModel, mcsModel,
                                   pixelGroupMap, CMInormal,
                                   gmcsGrowthLimits, gmcsMortLimits, gmcsMinAge){
  if (is.null(CMI) & is.null(ATA)) {
    message(paste("Missing climate data needed to run LandR.CS - consider running modules gmcsDataPrep and PSP_Clean",
                  "if you were expecting climate impacts for this year"))
    return(data.table('mortPred' = 100, 'growthPred' = 100))
  }

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

  climateMatch <- climateMatch[!is.na(pixelGroup)] #Not all pixelGroups are in pixelGroupMap, because cohortData is a subset
  #Take the median climate for each pixel group as some pixelgroups occur across multiple climate raster pixels
  climValues <- climateMatch[, .("CMI" = median(CMI, na.rm = TRUE),
                                 "ATA" = median(ATA, na.rm = TRUE),
                                 "CMInormal" = median(CMInormal, na.rm = TRUE)), by = "pixelGroup"]

  cohortData[, logAge := log(age)]
  setkey(cohortData, pixelGroup)
  setkey(climValues, pixelGroup)

  #Join cohort Data with climate data
  predData <- cohortData[climValues]

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

  #restrict predictions to those above min stand age
  climateEffect[age < gmcsMinAge, mortPred := 100]
  climateEffect[age < gmcsMinAge, growthPred := 100]

  climateEffect <- climateEffect[cohortData[, .(pixelGroup, speciesCode, age)], on = c('pixelGroup', 'speciesCode', 'age')]
  #this is to fix any pixelGroups that were dropped by the na.omit of climData due to NA climate values
  climateEffect[is.na(growthPred), c('growthPred', 'mortPred') := .(100, 100)]

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
#' @param x this is a param
#' @param y I don't know what this is
#' @param w I don't know what these are
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
