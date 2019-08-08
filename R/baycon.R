#' Bayesian deconvolution with Stan
#'
#' @export
#' @param numGenes numeric of number of genes
#' @param numCellTypes numeric of number of cell types
#' @param bulkExpression bulk expression matrix
#' @param sigMat signature matrix
#' @param mle a logical indicating whether we need mle or not
#' @param stanWarningCheck level at which Stan warnings must be handled, one of \code{'none', 'divergent', 'all'}
#' @param ... arguments to be passed to \code{rstan::sampling} (e.g. \code{chains, iter, init, verbose, refresh})
#' @return An object of type vector (list) with posterior mean estimates (and mle mean estimates)

baycon <- function(numGenes = nrow(p_bulkExpressionSimMat), numCellTypes = ncol(p_simSigMatTwo),
                   bulkExpression = p_bulkExpressionSimMat[,1:5], sigMat = p_simSigMatTwo,
                   mle = F, stanWarningCheck = "none", ...){

  mvnErrorDistList <- list()
  pEstimatesList <- list()
  sd2 <- list()

  for(i in 1:ncol(bulkExpression)) {

    # print(paste("******************************************************* ITERATION: ", i))
    standata <- list(numGenes = numGenes,numCellTypes = numCellTypes,
                     exprMixVec = bulkExpression[, i], sigMat = sigMat)

    nmfOut <- rstan::sampling(stanmodels$indSigmat, data = standata, ...);
    stanSumNmf <- as.data.frame(nmfOut)
    iters <- nrow(stanSumNmf)
    # do we want to restrict the number of iterations for which mean is computed below ?
    means <- apply(stanSumNmf, 2, mean)

    # store the purity point estimtes here (we're using the posterior mean. *Would mode be better??? Probably*)
    pEstimatesList[[i]] <- means[1:numCellTypes]

    if (mle) {
      # I need some error handling code because of general mlest buggyness.
      # Question: WHY does this just bomb sometimes ? Seems quite random. There must be a better approach here.
      tryCatch(
        # fit Multivariate normal by maximum liklihood, which will model the error.
        # QUESTION: How many samples do we pass to mlest
        mvnErrorDistList[[i]] <- mvnmle::mlest(data.matrix(as.data.frame(nmfOut)[1:iters, 1:numCellTypes])),
        error = function(cond) {
          message("******************************************************* ERROR HERE! ")
          message(cond)
          # Just literally try again until it works?
          i = i - 1
          # Seems to not work for no apparent reason sometimes.
          # Some stochastic component here?
          # Choose a return value in case of error
          return(NA)
        }
      )
    }

    # also capture the error estimates with a univariate normal distribution.
    # TODO?: we need to remove the warmup samples before estimating error
    # EDIT: this might be unnecessary as as.data.frame takes care of this
    # docs: Coerce the draws (without warmup) to an array, matrix or data frame. See as.array.stanfit.
    sd2[[i]] <- apply(stanSumNmf[, 1:numCellTypes], 2, sd)

    # print(paste("******************************************************* Done With Sample: ", i))
  }

  # save the estimates in a matrix
  propMat <- do.call(rbind, pEstimatesList)
  # QUESTION: Are we repeating this ?
  propMatStanEsts <- do.call(rbind, pEstimatesList)

  stan = list( mean = propMatStanEsts, errors = sd2)

  if (mle){
    # Pull out the proportion estimates that were quantified from the maximum likelihood estimates of the spatial data.
    # NB: Some of these dont work. I don't know why. Need to find a solution here.
    mvnPropList <- list()
    nullInd <- numeric()
    for(i in 1:ncol(bulkExpression))
    {
      mvnPropList[[i]] <- mvnErrorDistList[[i]]$muhat

      if(is.null(mvnErrorDistList[[i]]$muhat))
      {
        nullInd <- c(nullInd, i)
        message("NA detected")
        mvnPropList[[i]] <- rep(NA, numCellTypes)
      }
    }
    mvnPropEsts <- do.call(rbind, mvnPropList)
    mleEsts = list(mean = mvnPropEsts, mvnObj = mvnErrorDistList)
  }

  ret_value <- stan
  if (mle)
    ret_value <- list(stan = stan, mvn = mleEsts)

  return(ret_value)
}
