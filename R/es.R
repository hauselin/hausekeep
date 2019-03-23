#' @title Effect size conversion
#' @name es
#'
#' @description es is used to convert one effect size to other effect size measures. Currently supports Cohen's d, corelation r, r-squared, Cohen's f, odds ratio, log odds ratio, and area-under-curve (auc). Also available as a Shiny app here: http://escal.site
#'
#' @param d a numeric vector containing Cohen's d effect size(s)
#' @param r a numeric vector containing correlation r effect size(s)
#' @param R2 a numeric vector containing r-squared effect size(s)
#' @param f a numeric vector containing Cohen's f effect size(s)
#' @param oddsratio a numeric vector containing odds ratio effect size(s)
#' @param logoddsratio a numeric vector containing log odds ratio effect size(s)
#' @param auc a numeric vector containing area-under-curve effect size(s)
#' @param fishersz a numeric vector containing fisher's z effect size(s)
#' @param decimal a numeric vector indicating decimal places of output
#' @param msg a boolean indicating whether to show input effect size(s)
#'
#' @return A dataframe with converted effect sizes
#'
#' @details
#' Formulae for conversion
#' \cr \cr
#' \code{f = d / 2}
#' \cr
#' \code{r = d / sqrt(d^2 + 4)}
#' \cr
#' \code{d = (2 * r) / sqrt(1 - r^2)}
#' \cr
#' \code{R2 = r^2}
#' \cr
#' \code{oddsratio = exp(d / (sqrt(3) / pi))}
#' \cr
#' \code{logoddsratio = d / (sqrt(3) / pi)}
#' \cr
#' \code{auc = pnorm(d / sqrt(2), 0, 1)}
#' \cr
#' \code{fishers z = 0.5 * [log(1 + r) - log(1 - r)]}
#' \cr
#'
#' @note All conversions assume equal-sized groups. Effect size conventions:
#' \cr
#' Cohen's d: 0.20 (small), 0.50 (medium), .80 (large) (Cohen, 1992)
#' \cr
#' correlation r: .10 (small), .30 (medium), .50 (large)
#' \cr
#' R-squared: R2: .02 (small), .13 (medium), .26 (large)
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). Introduction to meta-analysis. Chichester, West Sussex, UK: Wiley.
#' \cr \cr
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.), Hillsdale, NJ: Erlbaum.
#' \cr \cr
#' Rosenthal, R. (1994). Parametric measures of effect size. In H. Cooper & L. V. Hedges (Eds.), The Handbook of Research Synthesis. New York, NY: Sage.
#' \cr \cr
#' Ruscio, J. (2008). A probability-based measure of effect size: Robustness to base rates and other factors. Psychological Methods, 13(1), 19-30. doi:10.1037/1082-989x.13.1.19
#'
#' @author Hause Lin
#'
#' @usage
#' es(d = NULL, r = NULL, R2 = NULL, f = NULL, oddsratio = NULL,
#' logoddsratio = NULL, auc = NULL, fishersz = NULL,
#' decimal = 3, msg = TRUE)
#'
#' @export
#' @examples
#' es(d = 0.3)
#' es(r = c(0.1, 0.3))
es <- function(d = NULL, r = NULL, R2 = NULL, f = NULL, oddsratio = NULL, logoddsratio = NULL, auc = NULL, fishersz = NULL, decimal = 3, msg = TRUE) {

  # effectsizes <- vector("list", 7) # list version
  effectsizes <- data.frame(matrix(NA, nrow = length(c(d, r, R2, f, oddsratio, logoddsratio, auc, fishersz)), ncol = 8)) # dataframe version
  names(effectsizes) <- c("d", "r", "R2", "f", "oddsratio", "logoddsratio", "auc", "fishersz")

  if (length(c(d, r, R2, f, oddsratio, logoddsratio, auc, fishersz)) < 1) {
    stop("Please specify one effect size!")
  }

  if (is.numeric(d)) {
    if (msg) {message(paste0("d: ", d, " ")) }
    effectsizes$d <- d
    effectsizes$f <- effectsizes$d / 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(r)) {
    if (msg) {message(paste0("r: ", r, " ")) }
    effectsizes$d <- (2 * r) / (sqrt(1 - r^2))
    effectsizes$r <- r
    effectsizes$f <- effectsizes$d / 2
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(f)) {
    if (msg) {message(paste0("f: ", f, " ")) }
    effectsizes$d <- f * 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$f <- f
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(R2)) {
    if (msg) {message(paste0("R2: ", R2, " ")) }
    effectsizes$r <- sqrt(R2)
    effectsizes$d <- (2 * effectsizes$r) / (sqrt(1 - effectsizes$r^2))
    effectsizes$f <- effectsizes$d / 2
    effectsizes$R2 <- R2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(oddsratio)) {
    if (msg) {message(paste0("odds ratio: ", oddsratio, " "))}
    effectsizes$d <- log(oddsratio) * (sqrt(3) / pi)
    effectsizes$f <- effectsizes$d / 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- oddsratio
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(logoddsratio)) {
    if (msg) {message(paste0("log odds ratio: ", logoddsratio, " ")) }
    effectsizes$logoddsratio <- logoddsratio
    effectsizes$d <- effectsizes$logoddsratio * (sqrt(3) / pi)
    effectsizes$f <- effectsizes$d / 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(auc)) { # also known as common language (CL) effect size statistic
    if (msg) {message(paste0("auc: ", auc, " ")) }
    effectsizes$auc <- auc
    effectsizes$d <- stats::qnorm(auc, 0, 1) * sqrt(2) # assumes equal sample size (Ruscio 2008)
    effectsizes$f <- effectsizes$d / 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(fishersz)) {
    if (msg) {message(paste0("fishersz: ", auc, " ")) }
    effectsizes$fishersz <- fishersz
    effectsizes$r <- (exp(effectsizes$fishersz / 0.5) - 1) / (exp(effectsizes$fishersz / 0.5) + 1)
    effectsizes$d <- (2 * effectsizes$r) / (sqrt(1 - effectsizes$r^2))
    effectsizes$f <- effectsizes$d / 2
    effectsizes$R2 <- R2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
  }

  round(effectsizes, decimal)
}
