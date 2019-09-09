QC of splatter simulation data
================
asif zubair
9/4/2019

## QC of `Splatter` simulations

Let’s either build a new simulation or load an old one

``` r
newSim = FALSE
if (newSim){
    source("make_splatter_sim.R")
} else {
  library(deconR)  
}
```

Take a look at some data:

``` r
if (newSim)
  SingleCellExperiment::rowData(splat_sim)
head(splat_sigMat)
```

    ##           Group1     Group2     Group3     Group4
    ## Gene1 35.3917197 47.9444444 47.9943503 46.6212121
    ## Gene2  0.3726115  0.4055556  0.4293785  0.4878788
    ## Gene3 14.5732484  8.7944444  9.5875706  8.9606061
    ## Gene4  6.2579618  8.2833333  9.2768362  8.8333333
    ## Gene5  1.0509554  1.4666667  1.4350282  1.4757576
    ## Gene6  0.8662420  1.3111111  1.1016949  1.2909091

and some plots:

``` r
if (newSim)
  scater::plotPCA(splat_sim, colour_by = "Group")
matplot(splat_sigMat, xlab = "Gene", ylab = "Expression")
```

![](splat_sim_qc_files/figure-gfm/plot-1.png)<!-- -->
