#' @title Compute z-score
#' @name zscore
#' @description Computes z-score (standardized scores)
#' @param x A vector
#' @param n use n or n-1 in the denominator when computing SD (default 0 for n; -1 for n = -1)
#' @return A vector
#'
#' @note Formula: (x - mean(x)) / SD(x)
#' @author Hause Lin
#' @export
#' @usage
#' zscore(x, n)
#' @examples
#' zscore(1:10)
#' zscore(1:10, n = 0) # default n = 0 (SD is computed using n)
#' zscore(1:10, n = -1) # n = -1 (SD is computed using n-1)
zscore <- function(x, n = 0) {
  x_without_na <- x[!is.na(x)]
  xmu <- mean(x_without_na)
  if (n == 0) {
    n <- length(x_without_na)
    stdev <- sqrt(sum((x_without_na - xmu)^2) / n)
    x_z <- (x_without_na - xmu) / stdev
  } else if (n == -1) {
    x_z <- (x_without_na - xmu) / sd(x_without_na)
  }
  return(x_z)
}
