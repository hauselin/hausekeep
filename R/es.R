#' Effect Size Conversion
#'
#' es is used to convert one effect size to other effect size measures
#'
#' Details to follow...
#'
#' @param d a numeric vector containing Cohen's d effect size(s)
#' @param r a numeric vector containing correlation r effect size(s)
#' @param R2 a numeric vector containing r-squared effect size(s)
#' @param f a numeric vector containing Cohen's f effect size(s)
#' @param oddsratio a numeric vector containing odds ratio effect size(s)
#' @param logoddsratio a numeric vector containing log odds ratio effect size(s)
#' @param auc a numeric vector containing area-under-curve effect size(s)
#' @param decimal a numeric vector indicating decimal places of output
#' @param msg a boolean indicating whether to show input effect size(s)
#'
#' @author Hause Lin
#'
#' @return A dataframe with difference effect sizes
#' @export
#'
#' @examples
#' es(d = 0.3)
#' es(r = c(0.1, 0.3))

es <- function(d = NULL, r = NULL, R2 = NULL, f = NULL, oddsratio = NULL, logoddsratio = NULL, auc = NULL, decimal = 3, msg = TRUE) {

    # effectsizes <- vector("list", 7) # list version
    effectsizes <- data.frame(matrix(NA, nrow = length(c(d, r, R2, f, oddsratio, logoddsratio, auc)), ncol = 7)) # dataframe version
    names(effectsizes) <- c("d", "r", "R2", "f", "oddsratio", "logoddsratio", "auc")
    # auc calculations might be off...

    if (length(c(d, r, R2, f, oddsratio, logoddsratio, auc)) < 1) {
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
        effectsizes$auc <- stats::pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(r)) {
        if (msg) {message(paste0("r: ", r, " ")) }
        effectsizes$d <- (2 * r) / (sqrt(1 - r^2))
        effectsizes$r <- r
        effectsizes$f <- effectsizes$d / 2
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- stats::pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(f)) {
        if (msg) {message(paste0("f: ", f, " ")) }
        effectsizes$d <- f * 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$f <- f
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- stats::pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(R2)) {
        if (msg) {message(paste0("R2: ", R2, " ")) }
        effectsizes$r <- sqrt(R2)
        effectsizes$d <- (2 * effectsizes$r) / (sqrt(1 - effectsizes$r^2))
        effectsizes$f <- effectsizes$d / 2
        effectsizes$R2 <- R2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- stats::pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(oddsratio)) {
        if (msg) {message(paste0("odds ratio: ", oddsratio, " "))}
        effectsizes$d <- log(oddsratio) * (sqrt(3) / pi)
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- oddsratio
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- stats::pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(logoddsratio)) {
        if (msg) {message(paste0("log odds ratio: ", logoddsratio, " ")) }
        effectsizes$logoddsratio <- logoddsratio
        effectsizes$d <- effectsizes$logoddsratio * (sqrt(3) / pi)
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$auc <- stats::pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(auc)) {
        if (msg) {message(paste0("auc: ", auc, " ")) }
        effectsizes$auc <- auc
        effectsizes$d <- stats::qnorm(auc, 0, 1)
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    }

    return(round(effectsizes, decimal))

}


# es(d = 0.3)
# es(d = 0.3, r = 0.2)
# es(d = c(0.2, 0.3, 0.4))
# es(d = c(0.2, 0.3))$r
# es(r = 0.5)
# es(f = 0.24)
# es(R2 = 0.6)
# es(oddsratio = 1.6)
# es(logoddsratio = 1.6)
# es(auc = .99)
# es()