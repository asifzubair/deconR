#'
#'
#'
#'
#'@return An object of list with the following names:
#'@return \code{Sensitivty, Specificty, MeasuredFalseDiscoveryRate,
#' numSigEQtls, FalseDiscoveriesIndex}

calcSensitivitySpecificityFdr <- function(estimatedBetas, pValues, realBetas, fdr = 0.2) {
  # Sensitivty (true positive rate) = number of true positives / number of real positives. (TP/P)
  # Calculate the number of "positives", i.e. the number of effects really different from zero.
  positives <- which(realBetas != 0)
  numPositives <- sum(realBetas != 0)

  # Calculate the number of True Positives,
  # i.e. the number of effects really different from zero that were identified by our model.
  correctDirectionality <- which((estimatedBetas * realBetas) > 0)
  isSignificant <- which(p.adjust(pValues, method="BH") < fdr)
  # things significant with correct directionality.
  numTruePositives <- sum((isSignificant %in% correctDirectionality) & (isSignificant %in% positives))

  sensitivity <- numTruePositives/numPositives

  # Specificity (true negative rate) = the number of true negatives / number of real negatives.
  # the number of real negatives
  negatives <- which(realBetas == 0)
  numNegatives <- sum(realBetas == 0)
  notSignificant <- which(p.adjust(pValues, method="BH") > fdr)
  # How many negatives were correctly assigned as "not significant"
  numTrueNegatives <- sum(notSignificant %in% negatives)
  specificity <- numTrueNegatives/numNegatives

  # Calculate the type I error rate, i.e. the type I error rate....!
  # i.e. the proportion of things called positive that are actually negative.
  MeasuredFalseDiscoveryRate <- sum(isSignificant %in% negatives) / length(isSignificant)
  FalseDiscoveriesIndex <- isSignificant[isSignificant %in% negatives]

  ret_value <- list(Sensitivty = sensitivity, Specificty = specificity,
                    MeasuredFalseDiscoveryRate = MeasuredFalseDiscoveryRate,
                    numSigEQtls = length(isSignificant), FalseDiscoveriesIndex = FalseDiscoveriesIndex)
  return(ret_value)
}
