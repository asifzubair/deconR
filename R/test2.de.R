#' Differential expression testing for two cell types
#'
#' This is a refactor of Paul's code.
#' Limited to two cell types.
#' NOT to be exported.
#'
#' @param y bulkExpression as a gene by samples matrix
#' @param genotype genotype labels for the samples
#' @param measProp2 estimated proportions of second cell type
#' @param sd2 standard deviation error of measured proportion
#' @param priorNormSD dunno what this is
#' @param numCelltypes number of cell types
#' @param n number of samples
#' @param log log progress
#' @param useHierarchichalModel use a hierarchical model instead
#' @param ... arguments to be passed to `rstan::sampling` (e.g. chains, iter (suggested > 8000), init, verbose, refresh)
#' @return list containing eQTL effect sizes for cancer and normal cell types


test2.de <- function(y, genotype, measProp2, sd2, priorNormSD = 100, numCellTypes, log = T, useHierarchichalModel = F, ...){

  # TODO: remove numCellTypes from function signature,
  # TODO: can compute it from one of the input matrices
  if (is.vector(y)) y <- t(as.data.frame(y))
  n <- ncol(y)

  ## Stan model, which accounts for differences in variance in underlying cell types,
  # and error associated with proportion estimates.
  # Variables.
  sumMcmcList <- list()
  betaCancer <- numeric()
  betaNormal <- numeric()
  pCancer <- numeric()
  pNormal <- numeric()

  # Compile the model outside the loop. Doing this inside the loop
  # (i.e. using "stan()" function causes segmentation faults (known bug))
  # model <- stan_model(file="error_model_2cellTypes.stan", model_name = "error_model_2cellTypes")

  # time this
  if (log) start.time <- Sys.time()
  for(i in 1:nrow(y))
  {
    if (log) print(paste("*********** ITERATION:", i, "*********** "))

    # NOTE: NOT SURE IF TRUE:
    # Note, the "noise" associated with the proportion measurements have not been measured in a reasonable way here,
    # i.e. the values used are just to test the implementation, but aren't relevant values in terms of the actual data.
    standata <- list(y = y[i,], genotype = genotype, measProp2 = measProp2, sd2 = sd2,
                      priorNormSD = 100, numCellTypes = numCellTypes, n = n)

    # run model on this data
    if (useHierarchichalModel){
      message("using hierarchical model")
      stanModOut <- sampling(object = stanmodels$errorModel2cellTypesHeirar, data = standata, ...)
    } else {
      message("using measurement error model")
      stanModOut <- sampling(object = stanmodels$errorModel2cellTypes, data = standata, ...)
    }

    # Summarize results.
    # There's prob a more efficient way of getting what we want out of this
    # TODO: Use rstan::get_posterior_mean function ?
    sumMcmcList[[i]] <- rstan::summary(stanModOut)$summary

    # recover the effect sizes for "cancer" and "normal"
    betaCancer[i] <- sumMcmcList[[i]]["beta1",1]
    betaNormal[i] <- sumMcmcList[[i]]["betaNormal",1]

    # Calculate some "p-values" (not technically p-values,
    # i.e. approximate the posterior distribution of the eQTL effects using a normal distribution)
    # The relationship between P-values and posterior probability is discussed in this paper,
    # basically, I think the approach below can be considered equivaent to a p-value assuming flat priors
    # (worth reading, although it mostly turns into a major critique of flat priors):
    # http://www.stat.columbia.edu/~gelman/research/published/pvalues3.pdf

    mcmcDraws <- rstan::extract(stanModOut)
    cancerDraws <- mcmcDraws[["beta1"]]
    normalDraws <- mcmcDraws[["betaNormal"]]

    cancerPUpperTail <- pnorm(0, mean(cancerDraws), sd(cancerDraws), lower.tail=F) #
    cancerPLowerTail <- pnorm(0, mean(cancerDraws), sd(cancerDraws), lower.tail=T)
    pCancer[i] <- (min(c(cancerPUpperTail, cancerPLowerTail)) * 2)

    normalPUpperTail <- pnorm(0, mean(normalDraws), sd(normalDraws), lower.tail=F) #
    normalPLowerTail <- pnorm(0, mean(normalDraws), sd(normalDraws), lower.tail=T)
    pNormal[i] <- (min(c(normalPUpperTail, normalPLowerTail)) * 2)

    if (log) print(paste("*********** Done with gene:", i, "*********** "))
  }

  if (log){
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  }

  ret_value <- list(stanObj = sumMcmcList, betaCancer = betaCancer, pCancer = pCancer, betaNormal = betaNormal, pNormal = pNormal)
  return(ret_value)

  }
