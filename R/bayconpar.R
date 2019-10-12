#' Bayesian deconvolution with Stan in parallel
#'
#' @export
#' @param bulkExpression bulk expression matrix
#' @param sigMat signature matrix
#' @param fit.mvn a logical indicating whether we want to fit a MVN to the posterior
#' @param useHyperPrior use a deconvolution model with hyper prior on the degrees of freedom
#' @param useOptim use the L-BFGS optimizer to determine posterior mode
#' @param useADVI use the variational bayes approach to estimate posterior mean and variances
#' @param ... arguments to be passed to \code{rstan::sampling} (e.g. \code{chains, iter, init, verbose, refresh})
#' @return An object of type list with posterior mean estimates (and MVN mean estimates),
#' alongwith associated error variances
bayconpar <- function(bulkExpression, sigMat, fit.mvn = F, useHyperPrior = F, useOptim = F, useADVI = F, retfunc = bayconResultClass, ...){

  if (is.numeric(bulkExpression)) bulkExpression <- as.data.frame(bulkExpression)
  if (any(useOptim, useADVI)) fit.mvn = FALSE
  stopifnot(nrow(bulkExpression) == nrow(sigMat))
  if (foreach::getDoParWorkers() == 1)
    message("parallel mode disabled")
  args = list(...)

  numGenes = nrow(bulkExpression)
  numCellTypes = ncol(sigMat)

  pEstimatesList <- list()
  var <- list()
  if (fit.mvn) mvnErrorDistList <- list()

  # message("parallel execution would follow metrics:",
  #        "\nName --> ", foreach::getDoParName(),
  #        "\nVersion --> ", foreach::getDoParVersion(),
  #        "\nWorkers --> ", foreach::getDoParWorkers(), "\n")

  if (useHyperPrior)
    model = stanmodels$indSigmatHyperprior
  else
    model = stanmodels$indSigmat

  parout <- foreach::foreach(i = 1:ncol(bulkExpression), .packages = "rstan") %dopar% {
    results <- retfunc()
    standata <- list(numGenes = numGenes, numCellTypes = numCellTypes,
                     exprMixVec = bulkExpression[, i], sigMat = sigMat)

    if (useOptim)
      nmfOut <- do.call(optimizing, c(list(object = model, data = standata, draws = 1000), args))
    else if (useADVI)
      nmfOut <- do.call(vb, c(list(object = model, data = standata, importance_resampling = T), args))
    else
      nmfOut <- do.call(sampling, c(list(object = model, data = standata), args))

    if (!(useOptim)){
      stanSumNmf <- as.data.frame(nmfOut)
      results$pEstimatesList <- colMeans(stanSumNmf[, 1:numCellTypes])
      results$var <- var(stanSumNmf[, 1:numCellTypes])
    } else {
      results$pEstimatesList <- nmfOut$par[1:numCellTypes]
      results$var <- apply(nmfOut$theta_tilde[, 1:numCellTypes], 2, sd)
    }

    if (fit.mvn)
      results$mvnErrorDistList <- mclust::mvn("XXX", stanSumNmf[, 1:numCellTypes], warn = T)
    return(results)
  }

  propMatStanEsts <- t(sapply(parout, "[[", "pEstimatesList"))
  var <- sapply(parout, "[[", "var", simplify = F)
  mvnErrorDistList <- sapply(parout, "[[", "mvnErrorDistList", simplify = F)

  # propMatStanEsts <- do.call(rbind, pEstimatesList)
  stanEsts = list(mean = propMatStanEsts, errors = var)

  if (fit.mvn){
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

#' Multi-result class for baycon
#'
#' Create class which holds multiple results for each loop iteration.
#' Each loop iteration populates three properties.
#' For a great tutorial on \code{S3} classes, see:
#' http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
#' @param pEstimatesList mean estimates from Stan
#' @param var variance estimates from Stan
#' @param mvnErrorDistList MVN approx to posterior
bayconResultClass <- function(pEstimatesList = NULL, var = NULL, mvnErrorDistList = NULL)
{
  me <- list(
    pEstimatesList = pEstimatesList,
    var = var,
    mvnErrorDistList = mvnErrorDistList
  )

  ## Set the name for the class
  class(me) <- append(class(me),"bayconResultClass")
  return(me)
}
