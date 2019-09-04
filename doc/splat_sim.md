Splatter Simulations
================
asif zubair
8/29/2019

``` r
library(deconR)
source("~/projects/decon/CIBERSORT/CIBERSORT.R")
options(mc.cores = parallel::detectCores())
```

## Data

    ## Num. of Cell Types:  4 
    ## Num. of Genes:  2000 
    ## Num. of Bulk Samples:  100

## Fit

``` r
t(splat_props[,1:5])
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,] 0.36 0.28 0.17 0.19
    ## [2,] 0.06 0.10 0.29 0.55
    ## [3,] 0.07 0.49 0.12 0.31
    ## [4,] 0.25 0.45 0.07 0.23
    ## [5,] 0.63 0.18 0.01 0.18

``` r
out <- baycon(bulkExpression = splat_bulkExpression, sigMat = splat_sigMat, useHyperPrior = T, refresh = 0, iter = 3000)
```

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#bulk-ess
    
    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#bulk-ess

    ## Warning: There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ## Warning: Examine the pairs() plot to diagnose sampling problems

    ## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
    ## Running the chains for more iterations may help. See
    ## http://mc-stan.org/misc/warnings.html#tail-ess

``` r
out$stan$mean[1:5,]
```

    ##      estimatedProportionsVecSimp[1] estimatedProportionsVecSimp[2]
    ## [1,]                     0.64269505                     0.06097352
    ## [2,]                     0.08859229                     0.29039826
    ## [3,]                     0.59208815                     0.09797526
    ## [4,]                     0.48535278                     0.16200457
    ## [5,]                     0.16650551                     0.31972726
    ##      estimatedProportionsVecSimp[3] estimatedProportionsVecSimp[4]
    ## [1,]                     0.24308049                     0.05325094
    ## [2,]                     0.12573515                     0.49527430
    ## [3,]                     0.21992679                     0.09000980
    ## [4,]                     0.12169113                     0.23095152
    ## [5,]                     0.02734933                     0.48641790

``` r
# let's also look at how linear regression does
pEstimate = list()
for (i in 1:ncol(splat_bulkExpression))
  pEstimate[[i]] <- coef(summary(lm(splat_bulkExpression[,i] ~ ., data = splat_sigMat)))[2:5,1]
lmEsts <- do.call(rbind, pEstimate)
lmEsts[1:5,]
```

    ##           Group1     Group2     Group3     Group4
    ## [1,]  0.70882958 -0.1775962  0.5226445 -0.0405499
    ## [2,] -0.09100362  0.4005454  0.2108778  0.4750463
    ## [3,]  0.27001641  0.2932814 -0.2384239  0.6968966
    ## [4,]  0.27194566  0.1866711 -0.3297387  0.8808299
    ## [5,]  0.07245644  0.4313767  0.3391111  0.1514394

``` r
cbs <- CIBERSORT(sig_matrix = splat_sigMat, mixture_file = splat_bulkExpression, perm  = 20)
cbs[1:5, ]
```

    ##          Group1    Group2    Group3    Group4 P-value Correlation
    ## Bulk1 0.3645666 0.1464047 0.3215129 0.1675158       0   0.9999103
    ## Bulk2 0.1457477 0.3365670 0.2362983 0.2813870       0   0.9999468
    ## Bulk3 0.1922081 0.2728235 0.1737676 0.3612008       0   0.9999541
    ## Bulk4 0.1406028 0.2001314 0.1922531 0.4670127       0   0.9999054
    ## Bulk5 0.1893763 0.2881135 0.2703636 0.2521466       0   0.9999606
    ##              RMSE
    ## Bulk1 0.013522494
    ## Bulk2 0.010337827
    ## Bulk3 0.009663857
    ## Bulk4 0.013857699
    ## Bulk5 0.008876614

![](splat_sim_files/figure-gfm/plot-1.png)<!-- -->

    ## Correlation: 0.09066219

![](splat_sim_files/figure-gfm/plot-2.png)<!-- -->

    ## Correlation: 0.06593122

![](splat_sim_files/figure-gfm/plot-3.png)<!-- -->

    ## Correlation: 0.07338582
