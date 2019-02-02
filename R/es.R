es <- function(d = NULL, r = NULL, R2 = NULL, f = NULL, oddsratio = NULL, logoddsratio = NULL, auc = NULL, decimal = 3, msg = TRUE) {
    
    # Last modified by Hause Lin 13-03-18 18:53 hauselin@gmail.com
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
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(r)) {
        if (msg) {message(paste0("r: ", r, " ")) }
        effectsizes$d <- (2 * r) / (sqrt(1 - r^2))
        effectsizes$r <- r
        effectsizes$f <- effectsizes$d / 2
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(f)) {
        if (msg) {message(paste0("f: ", f, " ")) }
        effectsizes$d <- f * 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$f <- f
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(R2)) {
        if (msg) {message(paste0("R2: ", R2, " ")) }
        effectsizes$r <- sqrt(R2)
        effectsizes$d <- (2 * effectsizes$r) / (sqrt(1 - effectsizes$r^2))
        effectsizes$f <- effectsizes$d / 2
        effectsizes$R2 <- R2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(oddsratio)) {
        if (msg) {message(paste0("odds ratio: ", oddsratio, " "))}
        effectsizes$d <- log(oddsratio) * (sqrt(3) / pi)
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- oddsratio
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(logoddsratio)) {
        if (msg) {message(paste0("log odds ratio: ", logoddsratio, " ")) }
        effectsizes$logoddsratio <- logoddsratio
        effectsizes$d <- effectsizes$logoddsratio * (sqrt(3) / pi)
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(auc)) {
        if (msg) {message(paste0("auc: ", auc, " ")) }
        effectsizes$auc <- auc
        effectsizes$d <- qnorm(auc, 0, 1)
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