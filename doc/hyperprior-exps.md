experiments with hyperprior
================
asif zubair
8/20/2019

## Motivation

We want to investigate the role/performance of the hyperprior when used
in the deconvolution algorithm. We will use some simple simulations to
do this.

## Normal distribution

Let’s say that the error is normally distributed.
Thus:

``` r
standata <- makeStanData(noisefunc = rnorm, n = 600, mean = 0, sd = .2)  
```

The output from the sampling algorithm with and without hyperprior is
below:

``` r
rstan::sampling(deconR:::stanmodels$indSigmat, data = standata, 
                verbose = F, refresh = 0);
```

    ## Inference for Stan model: indSigmat.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##                                  mean se_mean   sd   2.5%    25%    50%
    ## estimatedProportionsVecSimp[1]   0.75    0.00 0.02   0.71   0.73   0.75
    ## estimatedProportionsVecSimp[2]   0.25    0.00 0.02   0.22   0.24   0.25
    ## sigma                            0.13    0.00 0.00   0.12   0.12   0.13
    ## beta0                            0.00    0.00 0.01  -0.01  -0.01   0.00
    ## lp__                           845.85    0.03 1.24 842.73 845.30 846.15
    ##                                   75%  97.5% n_eff Rhat
    ## estimatedProportionsVecSimp[1]   0.76   0.78  3478    1
    ## estimatedProportionsVecSimp[2]   0.27   0.29  3478    1
    ## sigma                            0.13   0.14  3349    1
    ## beta0                            0.00   0.01  4668    1
    ## lp__                           846.75 847.29  2124    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Aug 29 15:10:19 2019.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

``` r
rstan::sampling(deconR:::stanmodels$indSigmatHyperprior, data = standata, 
                verbose = F, refresh = 0)
```

    ## Inference for Stan model: indSigmatHyperprior.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##                                  mean se_mean    sd   2.5%    25%    50%
    ## estimatedProportionsVecSimp[1]   0.75    0.00  0.02   0.72   0.74   0.75
    ## estimatedProportionsVecSimp[2]   0.25    0.00  0.02   0.21   0.24   0.25
    ## nu                              31.13    0.24 14.14  12.27  20.96  28.26
    ## sigma                            0.15    0.00  0.00   0.14   0.14   0.15
    ## beta0                            0.00    0.00  0.01  -0.01  -0.01   0.00
    ## lp__                           619.30    0.03  1.44 615.48 618.63 619.63
    ##                                   75%  97.5% n_eff Rhat
    ## estimatedProportionsVecSimp[1]   0.76   0.79  3716    1
    ## estimatedProportionsVecSimp[2]   0.26   0.28  3716    1
    ## nu                              38.13  65.99  3569    1
    ## sigma                            0.15   0.16  3471    1
    ## beta0                            0.00   0.01  4307    1
    ## lp__                           620.35 621.07  2096    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Aug 29 15:10:30 2019.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

We see that the `nu` value is close to 49.33. This is becasue of the the
normal error model that we use in the simulation. Essentially, with `df`
greater than 30, the t-distribution behaves like a normal. However, we
need to still be circumspect about how much influence the prior has on
the `nu` value.

Let’s do a few more tests to see if prior influence is indeed driving
this.

## t distribution with df = 4

``` r
standata <- makeStanData(noisefunc = deconR:::rgt, 
                         n = 600, df = 4, mu = 0, sigma = 1.2)
```

Let’s look at the output
now:

``` r
rstan::sampling(deconR:::stanmodels$indSigmat, data = standata, verbose = F, refresh = 0);
```

    ## Inference for Stan model: indSigmat.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##                                   mean se_mean   sd    2.5%     25%
    ## estimatedProportionsVecSimp[1]    0.82    0.00 0.11    0.58    0.75
    ## estimatedProportionsVecSimp[2]    0.18    0.00 0.11    0.01    0.09
    ## sigma                             1.02    0.00 0.04    0.94    0.99
    ## beta0                            -0.03    0.00 0.05   -0.13   -0.06
    ## lp__                           -430.88    0.03 1.27 -434.16 -431.48
    ##                                    50%     75%   97.5% n_eff Rhat
    ## estimatedProportionsVecSimp[1]    0.84    0.91    0.99  3910    1
    ## estimatedProportionsVecSimp[2]    0.16    0.25    0.42  3910    1
    ## sigma                             1.01    1.04    1.09  3359    1
    ## beta0                            -0.03    0.00    0.06  3338    1
    ## lp__                           -430.56 -429.95 -429.38  1701    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Aug 29 15:10:33 2019.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

``` r
rstan::sampling(deconR:::stanmodels$indSigmatHyperprior, data = standata, verbose = F, refresh = 0)
```

    ## Inference for Stan model: indSigmatHyperprior.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##                                   mean se_mean   sd    2.5%     25%
    ## estimatedProportionsVecSimp[1]    0.82    0.00 0.11    0.59    0.75
    ## estimatedProportionsVecSimp[2]    0.18    0.00 0.11    0.01    0.09
    ## nu                                4.41    0.02 0.90    3.03    3.80
    ## sigma                             1.03    0.00 0.06    0.92    0.99
    ## beta0                            -0.03    0.00 0.05   -0.13   -0.06
    ## lp__                           -674.34    0.04 1.51 -678.18 -675.09
    ##                                    50%     75%   97.5% n_eff Rhat
    ## estimatedProportionsVecSimp[1]    0.84    0.91    0.99  3385    1
    ## estimatedProportionsVecSimp[2]    0.16    0.25    0.41  3385    1
    ## nu                                4.27    4.88    6.58  1948    1
    ## sigma                             1.03    1.07    1.14  1953    1
    ## beta0                            -0.03    0.00    0.07  2677    1
    ## lp__                           -674.02 -673.21 -672.44  1611    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Aug 29 15:10:41 2019.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

We see that the hyperprior model recovers that `nu = 5.45` which is
pretty close to the actual value of 4.

## t distribution with df = 200

``` r
standata <- makeStanData(noisefunc = deconR:::rgt, 
                         n = 600, df = 200, mu = 0, sigma = 1.2)
```

Again, computing the output:

``` r
rstan::sampling(deconR:::stanmodels$indSigmat, data = standata, 
                verbose = F, refresh = 0);
```

    ## Inference for Stan model: indSigmat.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##                                   mean se_mean   sd    2.5%     25%
    ## estimatedProportionsVecSimp[1]    0.61    0.00 0.11    0.39    0.53
    ## estimatedProportionsVecSimp[2]    0.39    0.00 0.11    0.17    0.31
    ## sigma                             0.79    0.00 0.03    0.74    0.77
    ## beta0                             0.03    0.00 0.04   -0.05    0.00
    ## lp__                           -247.55    0.03 1.27 -250.96 -248.12
    ##                                    50%     75%   97.5% n_eff Rhat
    ## estimatedProportionsVecSimp[1]    0.61    0.69    0.83  3130    1
    ## estimatedProportionsVecSimp[2]    0.39    0.47    0.61  3130    1
    ## sigma                             0.79    0.81    0.85  3856    1
    ## beta0                             0.03    0.05    0.10  3520    1
    ## lp__                           -247.22 -246.61 -246.10  1969    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Aug 29 15:10:45 2019.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

``` r
rstan::sampling(deconR:::stanmodels$indSigmatHyperprior, data = standata, 
                verbose = F, refresh = 0)
```

    ## Inference for Stan model: indSigmatHyperprior.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##                                   mean se_mean    sd    2.5%     25%
    ## estimatedProportionsVecSimp[1]    0.61    0.00  0.11    0.39    0.53
    ## estimatedProportionsVecSimp[2]    0.39    0.00  0.11    0.17    0.32
    ## nu                               32.39    0.23 14.15   13.32   22.01
    ## sigma                             0.91    0.00  0.03    0.86    0.89
    ## beta0                             0.05    0.00  0.04   -0.02    0.03
    ## lp__                           -472.82    0.03  1.48 -476.47 -473.56
    ##                                    50%     75%   97.5% n_eff Rhat
    ## estimatedProportionsVecSimp[1]    0.61    0.68    0.83  3504    1
    ## estimatedProportionsVecSimp[2]    0.39    0.47    0.61  3504    1
    ## nu                               29.47   40.00   68.31  3668    1
    ## sigma                             0.91    0.94    0.98  3837    1
    ## beta0                             0.05    0.08    0.13  3912    1
    ## lp__                           -472.48 -471.72 -471.01  1819    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Aug 29 15:10:54 2019.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

Now we see tha the hyperprior model recovers that `nu = 25.56` which is
way off of the 200 true value.

## Conclusion

It seems that the prior model for `nu` that we use (`gamma(2, 0.1)`)
happens to be informative. This is why it struggles to recover really
high values of the degree of freedom. It is possible that the prior we
were using previously, `gamma(2, 0.01)`, was weakly informative and thus
we were able to recover both high and low values for the degrees of
freedom.
