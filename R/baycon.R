#' Bayesian deconvolution with Stan
#'
#' @export
#' @param numGenes numeric of number of genes
#' @param numCellTypes numeric of number of cell types
#' @param bulkExpression bulk expression matrix
#' @param sigMat signature matrix
#' @param mle a logical indicating whether we need mle or not
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains)
#' @return An object of type vector (list) with mean estimates from posterior (and mle estimates)

baycon <- function(numGenes = 1000, numCellTypes = 2, bulkExpression = p_bulkExpressionSimMat, sigmat = p_simSigMatTwo, mle = F, ...){

  mvnErrorDistList <- list()
  pEstimatesList <- list()
  sd1 <- numeric()
  sd2 <- numeric()

  data = list()

  for(i in 1:ncol(bulkExpression)) {

    standata <- list(numGenes = numGenes,
                     numCellTypes = numCellTypes,
                     exprMixVec = bulkExpression[, i],
                     sigMat = sigmat)

    nmfOut <- rstan::sampling(stanmodels$indSigmat, data = standata, chains = 3, iter = 2000, init = 0, ...);
    stanSumNmf <- summary(nmfOut)$summary

    # store the purity point estimtes here (we're using the posterior mean. *Would mode be better??? Probably*)
    pEstimatesList[[i]] <- stanSumNmf[c(1,2), 1]

    if (mle) {
      # I need some error handling code because of general mlest buggyness.
      # Question: WHY does thisjust bomb sometimes????? Seems quite random. There must be a better approach here.
      tryCatch(
        # fit Multivariate normal by maximum liklihood, which will model the error.
        mvnErrorDistList[[i]] <- mvnmle::mlest(data.matrix(as.data.frame(nmfOut)[1:numGenes, 1:numCellTypes])),
        error = function(cond) {
          message("***ERROR HERE***")
          message(cond)
          i = i - 1 # Just literally try again until it works??? Seems to not work for no apparent reason sometimes. Some stochastic component here?
          # Choose a return value in case of error
          return(NA)
        }
      )
    }

    # also capture the error estimates with a univariate normal distribution.
    sd1[i] <- sd(as.data.frame(nmfOut)[1:1000, 1])
    sd2[i] <- sd(as.data.frame(nmfOut)[1:1000, 2])
    print(paste("*********** ITERATION: ", i))
  }

  # save the estimates in a matrix
  propMat <- do.call(rbind, pEstimatesList)
  propMatStanEsts <- do.call(rbind, pEstimatesList)

  if (mle){
    # Pull out the proportion estimates that were quantified from the maximum liklihood estimtes of the spatial data.
    # NB: Some of these dont work. I don't know why. Need to find a solution here.
    mvnPropList <- list()
    nullInd <- numeric()
    for(i in 1:1000)
    {
      mvnPropList[[i]] <- mvnErrorDistList[[i]]$muhat

      if(is.null(mvnErrorDistList[[i]]$muhat))
      {
        nullInd <- c(nullInd, i)
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
