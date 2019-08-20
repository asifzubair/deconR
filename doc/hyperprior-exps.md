experiments with hyperprior
================
asif zubair
8/20/2019

## Motivation

We want to investigate the role/performance of the hyperprior when used
in the deconvolution algorithm. We will use some simple simulations to
do this.

## Normal distribution

Let’s say that the error is normally distributed. Thus:

``` r
  cancerExpressionMat[,i] <- seq(1, 1.599, .001) + rnorm(600, 0, .2)
  normalExpressionMat[,i] <- seq(1.599, 1, -.001) + rnorm(600, 0, .2)
  p_bulkExpressionSimMat[, i] <- (cancerExpressionMat[,i] * p_theProp[i]) + 
    (normalExpressionMat[,i] * p_propInv[i])
  p_simSigMatTwo <- cbind(p_cancerSig, p_normalSig)
  standata <- list(numGenes = numGenes, numCellTypes = numCellTypes, 
                   exprMixVec = p_bulkExpressionSimMat[,1], sigMat = p_simSigMatTwo)
```

The output from the sampling algorithm with and without hyperprior is
below:

``` r
> rstan::sampling(stanmodels$indSigmat, data = standata, verbose = F, refresh = 0);
Inference for Stan model: indSigmat.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                                 mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
estimatedProportionsVecSimp[1]   0.70    0.00 0.02   0.66   0.69   0.70   0.71   0.74  3247    1
estimatedProportionsVecSimp[2]   0.30    0.00 0.02   0.26   0.29   0.30   0.31   0.34  3247    1
sigma                            0.13    0.00 0.01   0.12   0.13   0.13   0.14   0.14  3219    1
beta0                            0.00    0.00 0.01  -0.01   0.00   0.00   0.01   0.01  4090    1
lp__                           815.33    0.03 1.29 811.99 814.76 815.66 816.25 816.79  1890    1

Samples were drawn using NUTS(diag_e) at Tue Aug 20 16:02:55 2019.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
> rstan::sampling(stanmodels$indSigmatHyperprior, data = standata, verbose = F, refresh = 0)
Inference for Stan model: indSigmatHyperprior.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                                 mean se_mean     sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
estimatedProportionsVecSimp[1]   0.70    0.00   0.02   0.66   0.69   0.70   0.71   0.74  3816    1
estimatedProportionsVecSimp[2]   0.30    0.00   0.02   0.26   0.29   0.30   0.31   0.34  3816    1
nu                             223.59    2.08 138.29  47.10 122.38 192.70 292.12 564.54  4421    1
sigma                            0.16    0.00   0.00   0.15   0.16   0.16   0.16   0.17  3933    1
beta0                            0.00    0.00   0.01  -0.01   0.00   0.00   0.01   0.02  4562    1
lp__                           597.13    0.03   1.36 593.64 596.46 597.46 598.12 598.82  2100    1

Samples were drawn using NUTS(diag_e) at Tue Aug 20 16:03:13 2019.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1)
```

We see that the `nu` value is close to 223. This shouldn’t surprise us
because of the normal error model that we use in the simulation.
However, we need to still be circumspect that the estimated `nu` value
could just be prior influence.

Let’s do a few more tests to see if prior influence is indeed driving
this.

## t distribution with df = 4

``` r
  cancerExpressionMat[,i] <- seq(1, 1.599, .001) + rt(600, 4)
  normalExpressionMat[,i] <- seq(1.599, 1, -.001) + rt(600, 4)
```

Let’s look at the output
now:

``` r
> rstan::sampling(stanmodels$indSigmat, data = standata, verbose = F, refresh = 0);
Inference for Stan model: indSigmat.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                                  mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
estimatedProportionsVecSimp[1]    0.65    0.00 0.10    0.44    0.58    0.65    0.72    0.85  3304    1
estimatedProportionsVecSimp[2]    0.35    0.00 0.10    0.15    0.28    0.35    0.42    0.56  3304    1
sigma                             0.74    0.00 0.03    0.69    0.72    0.74    0.76    0.80  3481    1
beta0                             0.02    0.00 0.04   -0.05    0.00    0.02    0.05    0.09  3682    1
lp__                           -239.91    0.03 1.27 -243.21 -240.50 -239.58 -238.97 -238.44  1710    1

Samples were drawn using NUTS(diag_e) at Tue Aug 20 16:06:25 2019.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
> rstan::sampling(stanmodels$indSigmatHyperprior, data = standata, verbose = F, refresh = 0)
Inference for Stan model: indSigmatHyperprior.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                                  mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
estimatedProportionsVecSimp[1]    0.65    0.00 0.10    0.44    0.58    0.65    0.72    0.85  3460    1
estimatedProportionsVecSimp[2]    0.35    0.00 0.10    0.15    0.28    0.35    0.42    0.56  3460    1
nu                                4.81    0.02 0.88    3.38    4.20    4.69    5.32    6.80  2721    1
sigma                             0.77    0.00 0.04    0.70    0.74    0.77    0.79    0.84  2739    1
beta0                             0.02    0.00 0.04   -0.05   -0.01    0.02    0.04    0.09  3417    1
lp__                           -482.65    0.04 1.51 -486.56 -483.37 -482.32 -481.55 -480.80  1534    1

Samples were drawn using NUTS(diag_e) at Tue Aug 20 16:06:37 2019.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```

We see that the hyperprior model recovers that `nu = 4.81` which is
pretty close to the actual value of 4.

## t distribution with df = 200

``` r
  cancerExpressionMat[,i] <- seq(1, 1.599, .001) + rt(600, 200)
  normalExpressionMat[,i] <- seq(1.599, 1, -.001) + rt(600, 200)
```

Again, computing the
output:

``` r
> rstan::sampling(stanmodels$indSigmat, data = standata, verbose = F, refresh = 0);
Inference for Stan model: indSigmat.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                                  mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
estimatedProportionsVecSimp[1]    0.69    0.00 0.09    0.50    0.62    0.69    0.75    0.86  3388    1
estimatedProportionsVecSimp[2]    0.31    0.00 0.09    0.14    0.25    0.31    0.38    0.50  3388    1
sigma                             0.66    0.00 0.02    0.61    0.64    0.65    0.67    0.70  3298    1
beta0                             0.03    0.00 0.03   -0.03    0.01    0.03    0.05    0.10  3756    1
lp__                           -135.96    0.03 1.25 -139.01 -136.55 -135.65 -135.04 -134.51  1974    1

Samples were drawn using NUTS(diag_e) at Tue Aug 20 16:09:49 2019.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
> rstan::sampling(stanmodels$indSigmatHyperprior, data = standata, verbose = F, refresh = 0)
Inference for Stan model: indSigmatHyperprior.
4 chains, each with iter=2000; warmup=1000; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=4000.

                                  mean se_mean     sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
estimatedProportionsVecSimp[1]    0.71    0.00   0.09    0.53    0.65    0.71    0.78    0.89  4109    1
estimatedProportionsVecSimp[2]    0.29    0.00   0.09    0.11    0.22    0.29    0.35    0.47  4109    1
nu                              212.73    2.16 137.62   39.46  112.63  181.53  280.95  557.26  4062    1
sigma                             0.78    0.00   0.02    0.74    0.77    0.78    0.80    0.83  4423    1
beta0                             0.03    0.00   0.03   -0.04    0.00    0.03    0.05    0.09  3910    1
lp__                           -356.62    0.03   1.43 -360.19 -357.33 -356.28 -355.56 -354.84  1765    1

Samples were drawn using NUTS(diag_e) at Tue Aug 20 16:10:01 2019.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```

Now we see tha the hyperprior model recovers that `nu = 212`.

## Conclusion

It seems that the hyperprior model learns very well that true
degrees-of-freedom for the error model. It is another question if this
helps to estimate the coefficients. As such it should, but in this
trivial example it doesn’t appear to be. It makes one think if we could
use a simple normal model instead of the t-distributed model for the
errors.
