#' Sample from a simplex
#'
#' This function samples from a simplex using a trick found here - http://blog.geomblog.org/2005/10/sampling-from-simplex.html
#' @param n number of dimensions
#' @param safe \code{logical}, if \code{TRUE} ensures there are no zeroes after rounding
#' @return a vector of RVs sampled from a unit simplex

sample.simplex <- function(n = 4, safe = T){
  out <- rep(0.00, n)
  while (any(out == 0.00)){
    rv.exp <- rexp(n)
    N <- sum(rv.exp)
    out <- round(rv.exp/N, 2)
    if (!safe)
      break
  }
  return(out)
}
