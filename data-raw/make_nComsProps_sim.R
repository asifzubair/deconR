## code to prepare `make_nComsProps_sim` dataset goes here

# load("sys.Rda")
set.seed(12345)

# exctract 01A (primary site) breast cancer samples.
nComsProps_brca <- nComsProps[nComsProps[,2] == "BRCA" & substring(nComsProps[,1], 14, 16) == "01A", "CPE"]
nComsProps_brca_noNa <- nComsProps_brca[!is.na(nComsProps_brca)]

# Proportions that I will use in subsequent analysis
# select 1000 samples at random (without replacement).
p_theProp <- sample(nComsProps_brca_noNa, 1000)
p_propInv <- (1-p_theProp)

# Create a simulated cancer datasets.
# 600 genes, 1000 patients, genes expression goes from 1:10 in even steps with gaussian noise.
cancerExpressionMat <- numeric(1000*600)
dim(cancerExpressionMat) <- c(600, 1000)
normalExpressionMat  <- numeric(1000*600)
dim(normalExpressionMat) <- c(600, 1000)
p_bulkExpressionSimMat <-numeric(1000*600)
dim(p_bulkExpressionSimMat) <- c(600, 1000)

for(i in 1:1000)
{
  cancerExpressionMat[,i] <- seq(1, 1.599, .001) + rnorm(600, 0, .2)
  normalExpressionMat[,i] <- seq(1.599, 1, -.001) + rnorm(600, 0, .2)
  p_bulkExpressionSimMat[, i] <- (cancerExpressionMat[,i] * p_theProp[i]) + (normalExpressionMat[,i] * p_propInv[i])
}

# the values before noise was added
p_cancerSig <- seq(1, 1.599, .001)
# QUESTION: why do we take the mean here ?
# Did you mean to do this:
# p_normalSig <- seq(1.599, 1, -.001) ?
p_normalSig <- apply(normalExpressionMat, 1, mean)
p_simSigMatTwo <- cbind(p_cancerSig, p_normalSig)

# p_ for paul_
usethis::use_data(p_theProp, p_propInv, p_cancerSig, p_normalSig, p_simSigMatTwo, p_bulkExpressionSimMat, overwrite = T)
