#' Function to make a simple simulated dataset for differential expression testing.
#'
#' The dataset created has 600 genes x 1000 samples and needs the true proportion values.
#' This works for two cell types - cancer and normal.
#' NOT to be exported.
#'
#' @param sdCancer std dev of cancer sample
#' @param sdNormal std dev of normal sample
#' @param seed set the seed for the simulation
#' @return An object of type list with following names:
#' @return \code{cancerExpressionSimMat, normalExpressionSimMat, bulkExpressionSimMat,
#' simulatedEffectSizesCancer, simulatedEffectSizesNormal}

makeSimDataset <- function(sdCancer = 1, sdNormal = 1, seed = 12345,
                           theProp = p_theProp, propInv = p_propInv)
{
  # Setting seed before running all of the below will ensure consistent result each time this is run
  set.seed(seed)

  # Create actual simulated Expression datasets of 1,000 samples and 600 genes....
  cancerExpressionSimMat <- numeric(600*1000)
  normalExpressionSimMat <- numeric(600*1000)
  bulkExpressionSimMat <- numeric(600*1000)
  normalEffectsIndependentOfCancer <- numeric(600*1000)
  dim(cancerExpressionSimMat) <- c(600, 1000)
  dim(normalExpressionSimMat) <- c(600, 1000)
  dim(bulkExpressionSimMat) <- c(600, 1000)
  dim(normalEffectsIndependentOfCancer) <- c(600, 1000)
  numSamps <- 250

  # Create 100 genes with a mixture of *different* eQTL in both Cancer and normal component.
  # Increment effect sizes from -1 to +1.
  # Create 100 simulated effect sizes....
  simulatedEffectSizesCancer <- seq(-.5, 0.49, .01)
  # randomize the above
  simulatedEffectSizesNormal <- sample(simulatedEffectSizesCancer)
  for(i in 1:100)
  {
    # NOTE: rnorm() will add noise and mean of 0 and sd of 1,
    # because the real data were standardized to a mean of 0 and sd of 1.
    cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdCancer),
                                    rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdCancer),
                                    rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdCancer))
    normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdNormal),
                                    rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdNormal),
                                    rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdNormal))
    # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
    bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv)
  }

  # Create 100 genes with an eQTL in cancer only
  simulatedEffectSizesCancer[101:200] <- seq(-.5, 0.49, .01)
  simulatedEffectSizesNormal[101:200] <- rep(0, 100)
  for(i in 101:200)
  {
    cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdCancer),
                                    rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdCancer),
                                    rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdCancer))
    normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdNormal),
                                    rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdNormal),
                                    rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdNormal))
    # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
    bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv)
  }

  # Create 100 genes with an eQTL in normal only
  simulatedEffectSizesCancer[201:300] <- rep(0, 100)
  simulatedEffectSizesNormal[201:300] <- seq(-.5, 0.49, .01)
  for(i in 201:300)
  {
    cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdCancer),
                                    rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdCancer),
                                    rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdCancer))
    normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdNormal),
                                    rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdNormal),
                                    rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdNormal))
    # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
    bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv)
  }

  # Create 100 genes with an eQTL in neither
  simulatedEffectSizesCancer[301:400] <- rep(0, 100)
  simulatedEffectSizesNormal[301:400] <- rep(0, 100)
  for(i in 301:400)
  {
    cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdCancer),
                                    rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdCancer),
                                    rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdCancer))
    normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdNormal),
                                    rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdNormal),
                                    rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdNormal))
    # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
    bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv)
  }

  # Create 100 genes with the SAME eQTL in Cancer and normal.
  simulatedEffectSizesCancer[401:500] <- seq(-.5, 0.49, .01)
  simulatedEffectSizesNormal[401:500] <- seq(-.5, 0.49, .01)
  for(i in 401:500)
  {
    cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdCancer),
                                    rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdCancer),
                                    rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdCancer))
    normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdNormal),
                                    rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdNormal),
                                    rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdNormal))
    # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
    bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv)
  }

  # Create 100 genes with SIMILAR eQTL in Cancer and Normal.
  # i.e. add some noise to the "cancer" eQTL effect size to create a "normal" effect size.
  simulatedEffectSizesCancer[501:600] <- seq(-.5, 0.49, .01)
  a <- seq(-.5, 0.49, .01) + rnorm(100, 0, .1)
  # Scale this so the effects won't be bigger than .5, to be consitent with everything else.
  aScaled <- (a / max(a)) *.5
  simulatedEffectSizesNormal[501:600] <- aScaled
  for(i in 501:600)
  {
    cancerExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdCancer),
                                    rep(1 + simulatedEffectSizesCancer[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdCancer),
                                    rep(1 + (simulatedEffectSizesCancer[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdCancer))
    normalExpressionSimMat[i,] <- c(rep(1, numSamps) + rnorm(numSamps, sd = sdNormal),
                                    rep(1 + simulatedEffectSizesNormal[i], (numSamps*2)) +
                                      rnorm(numSamps, sd = sdNormal),
                                    rep(1 + (simulatedEffectSizesNormal[i]*2), numSamps) +
                                      rnorm(numSamps, sd = sdNormal))
    # Combine the above to create a bulk expression data, assuming gene expression is additive based on the proportions.
    bulkExpressionSimMat[i, ] <- (cancerExpressionSimMat[i,] * theProp) + (normalExpressionSimMat[i,] * propInv)
  }

  ret_value <- list(cancerExpressionSimMat = cancerExpressionSimMat,
                    normalExpressionSimMat = normalExpressionSimMat,
                    bulkExpressionSimMat = bulkExpressionSimMat,
                    simulatedEffectSizesCancer = simulatedEffectSizesCancer,
                    simulatedEffectSizesNormal = simulatedEffectSizesNormal)
  return(ret_value)
  }
