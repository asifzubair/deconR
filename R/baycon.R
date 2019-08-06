#' Bayesian deconvolution with Stan
#'
#' @export
#' @param numGenes numeric of number of genes
#' @param numCellTypes numeric of number of cell types
#' @param bulkExpression bulk expression matrix
#' @param sigmat signature matrix
#' @param mle a logical indicating whether we need mle or not
#' @param ... arguments to be passed to `rstan::sampling` (e.g. chains, iter, init, verbose, refresh)
#' @return An object of type vector (list) with mean estimates from posterior (and mle estimates)

baycon <- function(numGenes = nrow(p_bulkExpressionSimMat), numCellTypes = ncol(p_simSigMatTwo),
                   bulkExpression = p_bulkExpressionSimMat, sigmat = p_simSigMatTwo,
                   mle = F, ...){

  mvnErrorDistList <- list()
  pEstimatesList <- list()

  for(i in 1:ncol(bulkExpression)) {

    print(paste("*********** ITERATION: ", i, " ***********"))
    standata <- list(numGenes = numGenes,
                     numCellTypes = numCellTypes,
                     exprMixVec = bulkExpression[, i],
                     sigMat = sigmat)

    nmfOut <- rstan::sampling(stanmodels$indSigmat, data = standata, ...);
    stanSumNmf <- as.data.frame(nmfOut)
    iters <- nrow(stanSumNmf)
    means <- apply(stanSumNmf, 2, mean)

    # store the purity point estimtes here (we're using the posterior mean. *Would mode be better??? Probably*)
    pEstimatesList[[i]] <- means[1:numCellTypes]

    if (mle) {
      # I need some error handling code because of general mlest buggyness.
      # Question: WHY does this just bomb sometimes ? Seems quite random. There must be a better approach here.
      tryCatch(
        # fit Multivariate normal by maximum liklihood, which will model the error.
        mvnErrorDistList[[i]] <- mvnmle::mlest(data.matrix(as.data.frame(nmfOut)[1:iters, 1:numCellTypes])),
        error = function(cond) {
          message("***ERROR HERE***")
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
    # TODO: we need to remove the warmup samples before estimating error
    # EDIT: this might be unnecessary as as.data.frame takes care of this
    # docs: Coerce the draws (without warmup) to an array, matrix or data frame. See as.array.stanfit.
    errors <- apply(stanSumNmf[1:iters, 1:numCellTypes], 2, sd)
    print(paste("*********** Done With Sample: ", i, " ***********"))
  }

  # save the estimates in a matrix
  propMat <- do.call(rbind, pEstimatesList)
  propMatStanEsts <- do.call(rbind, pEstimatesList)

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
        mvnPropList[[i]] <- c(NA, NA)
      }
    }
    mvnPropEsts <- do.call(rbind, mvnPropList)
  }

  ret_value = propMatStanEsts
  if (mle)
    ret_value = list(stan = propMatStanEsts, mle = mvnPropEsts)

  return(ret_value)
}
