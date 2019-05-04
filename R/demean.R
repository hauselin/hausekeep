#' @title De-mean or mean center variables
#' @name demean
#'
#' @description demean removes the mean or mean-centers a variable (excludes NA values)
#'
#' @param x a vector of numbers
#'
#' @return A vector of mean-centered numbers (mean = 0)
#'
#' @note This function is just a wrapper for x - mean(x, na.rm = T)
#'
#' @author Hause Lin
#'
#' @export
#'
#' @usage
#' demean(x)
#'
#' @examples
#' demean(1:3)
#'
demean <- function(x) {
  return(x - mean(x, na.rm = T))
}
