#' Curry function based on Stan debug level
#'
#' @param warn.level what level should the stan warnings be examined. one of \code{c('none', 'divergent', 'all')}
#' @return a Stan sampling function that handles warnings as indicated

