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


#' Samples from a generalized t-distribution
#'
#' Applies the specified location and scale paramters to the standart t-distribution.
#' The mean of this distribution is mu,
#' the variance is sigma^2*(df/(df - 2))
#' @param n number of samples to be generated
#' @param df degrees of freedom
#' @param mu the location parameter
#' @param sigma the scale parameter
#' @return a vector of generalized student-t distributed random variables

rgt <- function(n, df = 1, mu = 0, sigma = 1){
  out <- rt(n, df)
  out <- mu + sigma*out
  return(out)
}


#' Density of a generalized t-distribution
#'
#' This produces the density for a generalized t-distribution.
#' Details of the distribution are here -
#' https://en.wikipedia.org/wiki/Student\%27s_t-distribution#Generalized_Student's_t-distribution .
#' We want this density to match the one used in rstan here -
#' https://mc-stan.org/docs/2_18/functions-reference/student-t-distribution.html#probability-density-function-3
#'
#' @param x sample value
#' @param mu the location parameter
#' @param sigma the scale parameter
#' @return the density of the t-distribution at x

dgt <- function(x, df = 1, mu = 0, sigma = 1){
  stopifnot(df > 0)
  if (mu == 0 && sigma == 1) dt(x, df)

  out <- gamma((df + 1)/2)/gamma(df/2)
  out <- out * 1/sqrt(df*pi)
  out <- out * 1/sigma
  arg <- (1 + (1/df)*((x - mu)/sigma)**2)
  exponent <- -(df+1)/2
  out <- out*(arg**exponent)

  return(out)
}
