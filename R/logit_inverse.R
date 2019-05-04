#' @title Logistic (sigmoid or inverse logit) function
#' @name logit_inverse
#'
#' @description Applies inverse logit (sigmoid or logistic) transformation
#'
#' @param x a vector of numbers
#'
#' @return A vector ranging from 0 to 1
#'
#' @note Inverse logit can be expressed as exp(x) / (1 + exp(x)) or 1 / (1 + exp(-x))
#' \cr
#' See \url{https://hausetutorials.netlify.com/posts/2019-04-13-logistic-regression/#inverse-logit-and-logit-functions} for intro to logit and inverse logit functions.
#' @author Hause Lin
#'
#' @export
#'
#' @usage
#' logit_inverse(x)
#'
#' @examples
#' logit_inverse(-5:5)
#'
logit_inverse <- function(x) {
  return(1 / (1 + exp(-x))) # equivalent to exp(x) / (1 + exp(x))
}
