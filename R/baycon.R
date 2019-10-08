#' Bayesian deconvolution with Stan
#'
#' @export
#' @param bulkExpression bulk expression matrix
#' @param sigMat signature matrix
#' @param fit.mvn a logical indicating whether we want to fit a MVN to the posterior
#' @param useHyperPrior use a deconvolution model with hyper prior on the degrees of freedom
#' @param useOptim use the L-BFGS optimizer to determine posterior mode
#' @param useADVI use the variational bayes approach to estimate posterior mean and variances
#' @param stanWarningCheck level at which Stan warnings must be handled, one of \code{c('none', 'divergent', 'all')}
#' @param ... arguments to be passed to \code{rstan::sampling} (e.g. \code{chains, iter, init, verbose, refresh})
#' @return An object of type list with posterior mean estimates (and MVN mean estimates),
#' alongwith associated error variances
baycon <- function(bulkExpression, sigMat, fit.mvn = F, useHyperPrior = F, useOptim = F, useADVI = F, stanWarningCheck = 'none', ...){

  # TODO: this following sanity check is very basic,
  # TODO: need to hand situations in which either row or column vectors are passed
  if (is.numeric(bulkExpression)) bulkExpression <- as.data.frame(bulkExpression)
  if (any(useOptim, useADVI)) fit.mvn = FALSE
  stopifnot(nrow(bulkExpression) == nrow(sigMat))
  if (foreach::getDoParWorkers() == 1)
    message("parallel mode disabled")

  numGenes = nrow(bulkExpression)
  numCellTypes = ncol(sigMat)

  pEstimatesList <- list()
  var <- list()
  if (fit.mvn) mvnErrorDistList <- list()

  #TODO: ideally, I want to do some diagnostic checks based on Stan warnings.
  #TODO: will probably write a curried function that does diagnostic checks.

  if (useHyperPrior)
    model = stanmodels$indSigmatHyperprior
  else
    model = stanmodels$indSigmat

  for(i in 1:ncol(bulkExpression)) {
    # print(paste("******************************************************* ITERATION: ", i))
    standata <- list(numGenes = numGenes, numCellTypes = numCellTypes,
                     exprMixVec = bulkExpression[, i], sigMat = sigMat)

    if (useOptim)
      nmfOut <- rstan::optimizing(model, data = standata, draws = 1000, ...)
    else if (useADVI)
      nmfOut <- rstan::vb(model, data = standata, importance_resampling = T, ...)
    else
      nmfOut <- rstan::sampling(model, data = standata, ...)

    if (!(useOptim)){
      stanSumNmf <- as.data.frame(nmfOut)
      # store the purity point estimtes here (we're using the posterior mean. *Would mode be better??? Probably*)
      # QUESTION: do we want to restrict the number of iterations for which mean is computed below ?
      pEstimatesList[[i]] <- colMeans(stanSumNmf[, 1:numCellTypes])
      # also capture the error estimates with a univariate normal distribution.
      # TODO?: we need to remove the warmup samples before estimating error
      # EDIT: this might be unnecessary as as.data.frame takes care of this
      # docs: Coerce the draws (without warmup) to an array, matrix or data frame. See as.array.stanfit.
      # NOTE: this gives unbiased estimate of the covariance matrix
      var[[i]] <- var(stanSumNmf[, 1:numCellTypes])
    } else {
      pEstimatesList[[i]] <- nmfOut$par[1:numCellTypes]
      var[[i]] <- apply(nmfOut$theta_tilde[, 1:numCellTypes], 2, sd)
    }

    if (fit.mvn)
      mvnErrorDistList[[i]] <- mclust::mvn("XXX", stanSumNmf[, 1:numCellTypes], warn = T)
    # print(paste("******************************************************* Done With Sample: ", i))
  }

  # save the estimates in a matrix
  # QUESTION: Are we repeating this ? Commenting out below line for now
  # propMat <- do.call(rbind, pEstimatesList)
  propMatStanEsts <- do.call(rbind, pEstimatesList)
  stanEsts = list(mean = propMatStanEsts, errors = var)

  if (fit.mvn){
    # Pull out the proportion estimates that were quantified from the maximum likelihood estimates of the spatial data.
    mvnPropEsts <- do.call(cbind,
                           sapply(mvnErrorDistList, "[[", "parameters")["mean", ])
    mvnPropEsts <- t(mvnPropEsts)
    mvnEsts = list(mean = mvnPropEsts, mvnObj = mvnErrorDistList)
  }

  if (fit.mvn)
    ret_value <- list(stan = stanEsts, mvn = mvnEsts)
  else
    ret_value <- list(stan = stanEsts, mvn = NULL)

  return(ret_value)
  }
