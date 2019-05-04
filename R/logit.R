#' @title Logit function
#' @name logit
#'
#' @description Applies logit transformation
#'
#' @param p a vector of probability values (0 to 1)
#'
#' @return A vector
#'
#' @note Formula: log(p / (1 - p))
#' \cr Returns NA values if input values are not between 0 and 1.
#' \cr
#' See \url{https://hausetutorials.netlify.com/posts/2019-04-13-logistic-regression/#inverse-logit-and-logit-functions} for intro to logit and inverse logit functions.
#'
#' @author Hause Lin
#'
#' @export
#'
#' @usage
#' logit(p)
#'
#' @examples
#' logit(seq(-2, 2, length.out = 100))
#'
logit <- function(p) {
  return(log(p / (1 - p))) # equivalent to exp(x) / (1 + exp(x))
}
