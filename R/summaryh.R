#' @title Generate formatted results and effect sizes for manuscripts
#' @name summaryh
#'
#' @description After fitting models, use summaryh() in place of summary() to get APA (American Psychological Association) formatted output that also includes effect size estimates for each effect (r effect size). Currently supports lm, glm, aov, anova, lmer, lme, t-test, chisq.test, and cor.test. Unfortunately, this function won't write your entire results section for you (yet).
#'
#' @param model a fitted model
#' @param decimal round output to decimal places
#' @param showTable show results in table format (returns list)
#' @param showEffectSizesTable show other effect sizes computed using es function
#' @param ... further arguments passed to or from other methods
#' @return A datatable or a list with datatables (if showTable = TRUE or showEffectSizesTable = TRUE)
#'
#' @import data.table
#' @note
#' Cohen's d: 0.20 (small), 0.50 (medium), .80 (large) (Cohen, 1992)
#' \cr
#' correlation r: .10 (small), .30 (medium), .50 (large)
#' \cr
#' R-squared: R2: .02 (small), .13 (medium), .26 (large)
#' @author Hause Lin
#' @export
#'
#' @usage
#' summaryh(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...)
#'
#' @examples
#' summaryh(lm(mpg ~ qsec, mtcars))
#' summaryh(aov(mpg ~ gear, mtcars))
#' summaryh(cor.test(mtcars$mpg, mtcars$gear), showEffectSizesTable = TRUE)
#' summaryh(t.test(mpg ~ vs, mtcars), showTable = TRUE)
#' summaryh(glm(vs ~ 1, mtcars, family = "binomial"), showTable = TRUE)
summaryh <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {
    UseMethod("summaryh")
}

summaryh.default <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {
  return(summary(model))
}

#### summaryh class methods
# anova
#' @importFrom sjstats cohens_f
#' @export
summaryh.aov <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }

    # example output: F(3, 10) = 39, p < .001, r = 0.32
    if (class(model)[1] == "anova") {
        estimates <- data.frame(model) # get estimates and put in dataframe
    } else if (class(model)[1] == "aov") {
        estimates <- data.frame(stats::anova(model)) # get estimates and put in dataframe
    }
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    estimates$dfResid <- estimates["Residuals", "df"]  # get model degrees of freedom
    colnames(estimates) <- c('df', 'sum.sq', 'mean.sq', 'f.value', 'p.value', 'df.resid') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL
    estimates <- estimates[estimates$term != "Residuals", ]

    # effect sizes
    esCohensf <- cohens_f(model) # calculate Cohen's f
    if (is.data.frame(esCohensf)) { # sometimes output is a dataframe; if so, extract cohens.f variable
        estimates$es.f <- esCohensf$cohens.f
    } else {
        estimates$es.f <- esCohensf
    }

    estimates$es.r <- es(f = estimates$es.f, msg = F, decimal = decimal)$r

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

    # leave df as integers
    estimatesRound$df <- round(estimates$df)
    estimatesRound$df.resid <- round(estimates$df.resid)

    formattedOutput <- paste0("F(", estimatesRound$df, ", ", estimatesRound$df.resid, ")",
                              " = ", estimatesRound$f.value,
                              ", p ", estimatesRound$p.value,
                              ", r = ", estimatesRound$es.r)

    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
    }

    formattedOutputDf <- data.table(term = as.character(estimates$term),
                                    results = as.character(formattedOutput))

    outputList <- list(results = formattedOutputDf)

    if (showTable) {

        # format table nicely
        estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
        estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)
        outputList$results2 <- estimatesOutput
    }

    if (showEffectSizesTable) {

        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal, msg = F)
        outputList$effectSizes <- data.table(term = as.character(estimates$term), effectSizes)

    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}

# anova
#' @export
summaryh.anova <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {

  # ensure significant digits with sprintf
  digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
  if (decimal <= 2) {
    pdigits <- paste0("%.", 3, "f")
  } else {
    pdigits <- paste0("%.", decimal, "f")
  }

  # example output: F(3, 10) = 39, p < .001, r = 0.32
  if (class(model)[1] == "anova") {
    estimates <- data.frame(model) # get estimates and put in dataframe
  } else if (class(model)[1] == "aov") {
    estimates <- data.frame(stats::anova(model)) # get estimates and put in dataframe
  }
  effectNames <- rownames(estimates) # get names of effects
  colnames(estimates) <- tolower(colnames(estimates))

  if (names(model)[1] == 'Sum Sq') { # if output format resembles anova table from a lmer model...
    estimates['dfResid'] <- estimates["dendf"]  # get model degrees of freedom
    colnames(estimates) <- c('sum.sq', 'mean.sq', 'df', 'df.resid', 'f.value', 'p.value', 'df.resid') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL
    estimates <- estimates[estimates$term != "Residuals", ]
    esCohensf <- 0
  } else {
    estimates$dfResid <- estimates["Residuals", "df"]  # get model degrees of freedom
    colnames(estimates) <- c('df', 'sum.sq', 'mean.sq', 'f.value', 'p.value', 'df.resid') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL
    estimates <- estimates[estimates$term != "Residuals", ]
    # effect size
    esCohensf <- cohens_f(model) # calculate Cohen's f
  }

  # effect sizes
  if (is.data.frame(esCohensf)) { # sometimes output is a dataframe; if so, extract cohens.f variable
    estimates$es.f <- esCohensf$cohens.f
  } else {
    estimates$es.f <- esCohensf
  }

  estimates$es.r <- es(f = estimates$es.f, msg = F, decimal = decimal)$r

  # make a copy of estimates and convert to correct dp
  estimatesCopy <- estimates[, -1]
  estimatesRound <- estimatesCopy
  estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
  estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
  estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
  estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
  # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

  # fix p values
  estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
  estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
  estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

  # leave df as integers
  estimatesRound$df <- round(estimates$df)
  estimatesRound$df.resid <- round(estimates$df.resid)

  formattedOutput <- paste0("F(", estimatesRound$df, ", ", estimatesRound$df.resid, ")",
                            " = ", estimatesRound$f.value,
                            ", p ", estimatesRound$p.value,
                            ", r = ", estimatesRound$es.r)

  # convert hyphens to minus (only possible on UNIX systems)
  if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
    formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
  }

  formattedOutputDf <- data.table(term = as.character(estimates$term),
                                  results = as.character(formattedOutput))

  outputList <- list(results = formattedOutputDf)

  if (showTable) {

    # format table nicely
    estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
    estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)
    outputList$results2 <- estimatesOutput
  }

  if (showEffectSizesTable) {

    # get all other effect sizes
    effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal, msg = F)
    outputList$effectSizes <- data.table(term = as.character(estimates$term), effectSizes)

  }
  options(scipen = 0) # enable scientific notation
  if (length(outputList) > 1) {
    return(outputList)
  } else {
    return(formattedOutputDf)
  }

}

# lm
#' @export
summaryh.lm <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, confInterval = NULL, ...) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }

    # example output: b = −2.88, SE = 0.32, t(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(stats::coef(summary(model))) # get estimates and put in dataframe
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    estimates$df <- stats::df.residual(model) # get model degrees of freedom
    estimates <- estimates[, c(1, 2, 5, 3, 4)] # sort columns
    colnames(estimates) <- c('estimate', 'std.error', 'df', 'statistic', 'p.value') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL

    if (!is.null(confInterval)) {
        confIntervalChar <- paste0(confInterval * 100, "%")
        confIntervals <- stats::confint(model, level = confInterval)
        if (nrow(estimates) == 1) {
            names(confIntervals) <- c('ciLower', 'ciUpper')
            confIntervals <- data.frame(ciLower = confIntervals[1], ciUpper = confIntervals[2])
            if (nrow(confIntervals) > 1) {
                confIntervals <- confIntervals[rownames(confIntervals) %in% estimates$term, ]
            }
        } else {
            colnames(confIntervals) <- c('ciLower', 'ciUpper')
            confIntervals <- confIntervals[rownames(confIntervals) %in% estimates$term, ]
        }
        rownames(confIntervals) <- NULL
        estimates <- cbind(estimates, confIntervals)
    }

    # effect sizes
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d

    # R2
    estimates$es.r.squared <- summary(model)$r.squared
    estimates$es.adj.r.squared <- summary(model)$adj.r.squared

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    if (!is.na(sum(estimatesRound$statistic))) {
      estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
      estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
      estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
      estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
      # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

      # fix p values
      estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
      estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
      estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
    }

    # leave df as integers
    estimatesRound$df <- round(estimates$df)

    if (!is.null(confInterval)) { # report CIs, not SE
        formattedOutput <- paste0("b = ", estimatesRound$estimate,
                                  ", ", confIntervalChar, " CI [", estimatesRound$ciLower, " ", estimatesRound$ciUpper, "]",
                                  ", t(", estimatesRound$df, ")",
                                  " = ", estimatesRound$statistic,
                                  ", p ", estimatesRound$p.value,
                                  ", r = ", estimatesRound$es.r)
    } else { # report SE, not CIs
        formattedOutput <- paste0("b = ", estimatesRound$estimate,
                                  ", SE = ", estimatesRound$std.error,
                                  ", t(", estimatesRound$df, ")",
                                  " = ", estimatesRound$statistic,
                                  ", p ", estimatesRound$p.value,
                                  ", r = ", estimatesRound$es.r)
    }

    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
    }

    formattedOutputDf <- data.table(term = as.character(estimates$term),
                                    results = as.character(formattedOutput))

    outputList <- list(results = formattedOutputDf)

    if (showTable) {

        # format table nicely
        estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
        estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)
        outputList$results2 <- estimatesOutput
    }

    if (showEffectSizesTable) {

        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal, msg = F)
        outputList$effectSizes <- data.table(term = as.character(estimates$term), effectSizes)

    }

    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}

# glm
#' @export
summaryh.glm <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, confInterval = NULL, ...) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }

    # example output: b = −2.88, SE = 0.32, z(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(stats::coef(summary(model))) # get estimates and put in dataframe
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    estimates$df <- stats::df.residual(model) # get model degrees of freedom
    colnames(estimates) <- c('estimate', 'std.error', 'statistic', 'p.value', 'df') # rename columns
    estimates <- data.frame(term = effectNames, estimates)
    rownames(estimates) <- NULL

    if (!is.null(confInterval)) {
        confIntervalChar <- paste0(confInterval * 100, "%")
        confIntervals <- stats::confint(model, level = confInterval)
        if (nrow(estimates) == 1) {
            names(confIntervals) <- c('ciLower', 'ciUpper')
            confIntervals <- data.frame(ciLower = confIntervals[1], ciUpper = confIntervals[2])
            if (nrow(confIntervals) > 1) {
                confIntervals <- confIntervals[rownames(confIntervals) %in% estimates$term, ]
            }
        } else {
            colnames(confIntervals) <- c('ciLower', 'ciUpper')
            confIntervals <- confIntervals[rownames(confIntervals) %in% estimates$term, ]
        }
        rownames(confIntervals) <- NULL
        estimates <- cbind(estimates, confIntervals)
    }

    # effect sizes
    estimates$es.oddsratio <- exp(estimates$estimate)
    estimates$es.r <- es(oddsratio = estimates$es.oddsratio, msg = F, decimal = decimal)$r # r
    estimates$es.d <- es(oddsratio = estimates$es.oddsratio, msg = F, decimal = decimal)$d # d

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

    # leave df as integers
    estimatesRound$df <- round(estimates$df)

    if (!is.null(confInterval)) { # report CIs, not SE
        formattedOutput <- paste0("b = ", estimatesRound$estimate,
                                  ", ", confIntervalChar, " CI [", estimatesRound$ciLower, " ", estimatesRound$ciUpper, "]",
                                  ", z(", estimatesRound$df, ")",
                                  " = ", estimatesRound$statistic,
                                  ", p ", estimatesRound$p.value,
                                  ", r = ", estimatesRound$es.r)
    } else { # report SE, not CIs
        formattedOutput <- paste0("b = ", estimatesRound$estimate,
                                  ", SE = ", estimatesRound$std.error,
                                  ", z(", estimatesRound$df, ")",
                                  " = ", estimatesRound$statistic,
                                  ", p ", estimatesRound$p.value,
                                  ", r = ", estimatesRound$es.r)
    }

    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
    }

    formattedOutputDf <- data.table(term = as.character(estimates$term),
                                    results = as.character(formattedOutput))

    outputList <- list(results = formattedOutputDf)

    if (showTable) {

        # format table nicely
        estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
        estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)
        outputList$results2 <- estimatesOutput
    }

    if (showEffectSizesTable) {

        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal)
        outputList$effectSizes <- data.table(term = as.character(estimates$term), effectSizes)

    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}

# glmer
#' @export
summaryh.glmerMod <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, confInterval = NULL, ...) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }

    # example output: b = −2.88, SE = 0.32, z(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(stats::coef(summary(model))) # get estimates and put in dataframe
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    estimates$df <- stats::df.residual(model) # get model degrees of freedom
    colnames(estimates) <- c('estimate', 'std.error', 'statistic', 'p.value', 'df') # rename columns
    estimates <- data.frame(term = effectNames, estimates)
    rownames(estimates) <- NULL

    if (!is.null(confInterval)) {
        confIntervalChar <- paste0(confInterval * 100, "%")
        confIntervals <- stats::confint(model, level = confInterval)
        if (nrow(estimates) == 1) {
            names(confIntervals) <- c('ciLower', 'ciUpper')
            confIntervals <- data.frame(ciLower = confIntervals[1], ciUpper = confIntervals[2])
            if (nrow(confIntervals) > 1) {
                confIntervals <- confIntervals[rownames(confIntervals) %in% estimates$term, ]
            }
        } else {
            colnames(confIntervals) <- c('ciLower', 'ciUpper')
            confIntervals <- confIntervals[rownames(confIntervals) %in% estimates$term, ]
        }
        rownames(confIntervals) <- NULL
        estimates <- cbind(estimates, confIntervals)
    }

    # effect sizes
    estimates$es.oddsratio <- exp(estimates$estimate)
    estimates$es.r <- es(oddsratio = estimates$es.oddsratio, msg = F, decimal = decimal)$r # r
    estimates$es.d <- es(oddsratio = estimates$es.oddsratio, msg = F, decimal = decimal)$d # d

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

    # leave df as integers
    estimatesRound$df <- round(estimates$df)

    if (!is.null(confInterval)) { # report CIs, not SE
        formattedOutput <- paste0("b = ", estimatesRound$estimate,
                                  ", ", confIntervalChar, " CI [", estimatesRound$ciLower, " ", estimatesRound$ciUpper, "]",
                                  ", z(", estimatesRound$df, ")",
                                  " = ", estimatesRound$statistic,
                                  ", p ", estimatesRound$p.value,
                                  ", r = ", estimatesRound$es.r)
    } else { # report SE, not CIs
        formattedOutput <- paste0("b = ", estimatesRound$estimate,
                                  ", SE = ", estimatesRound$std.error,
                                  ", z(", estimatesRound$df, ")",
                                  " = ", estimatesRound$statistic,
                                  ", p ", estimatesRound$p.value,
                                  ", r = ", estimatesRound$es.r)
    }

    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
    }

    formattedOutputDf <- data.table(term = as.character(estimates$term),
                                    results = as.character(formattedOutput))

    outputList <- list(results = formattedOutputDf)

    if (showTable) {

        # format table nicely
        estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
        estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)
        outputList$results2 <- estimatesOutput
    }

    if (showEffectSizesTable) {

        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal)
        outputList$effectSizes <- data.table(term = as.character(estimates$term), effectSizes)

    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}

# lmer without lmerTest
#' @export
summaryh.lmerMod <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {
    message("Please install/load lmerTest package and then refit your model!")
}

# lmer with lmerTest
#' @export
summaryh.merModLmerTest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }

    # example output: b = −2.88, SE = 0.32, t(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(stats::coef(summary(model))) # get estimates and put in dataframe
    if (ncol(estimates) < 5 & (class(model)[1] == "merModLmerTest")) {
        return(message("lmerTest failed to compute p values; use summary() to check; try refitting model with lme() from nlme package"))
    }
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    colnames(estimates) <- c('estimate', 'std.error', 'df', 'statistic', 'p.value') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL

    # effect size r (Kashdan & Steger, 2006)
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

    # leave df as integers
    estimatesRound$df <- round(estimates$df)

    formattedOutput <- paste0("b = ", estimatesRound$estimate,
                              ", SE = ", estimatesRound$std.error,
                              ", t(", estimatesRound$df, ")",
                              " = ", estimatesRound$statistic,
                              ", p ", estimatesRound$p.value,
                              ", r = ", estimatesRound$es.r)

    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
    }

    formattedOutputDf <- data.table(term = as.character(estimates$term),
                                    results = as.character(formattedOutput))

    outputList <- list(results = formattedOutputDf)

    if (showTable) {

        # effect size semi partial R (Edwards et al., 2008)
        anovaModel <- data.frame(stats::anova(model))
        colnames(anovaModel) <- tolower(colnames(anovaModel))
        Fs <- anovaModel$f.value # F-values for each effect (marginal = type 3 SS with Satterthwaite (requires lmerTest package))
        numDF <- anovaModel$numdf #numerator DFs
        denDF <- anovaModel$dendf #denominator DFs
        semiPartialREffect <- (numDF / denDF * Fs) / (1 + (numDF / denDF * Fs)) # effect sizes
        semiPartialREffect <- data.frame(term = rownames(anovaModel), es.partR2 = semiPartialREffect, stringsAsFactors = F)
        for (i in 1:nrow(semiPartialREffect)) {
            estimates$es.partR2[grepl(semiPartialREffect$term[i], estimates$term)] <- semiPartialREffect[i, "es.partR2"]
        }

        # piecewiseSEM (Nakagawa & Schielzeth, 2013)
        # rsquareds <- sem.model.fits(model)
        # estimates$es.R2marginal <- rsquareds$Marginal
        # estimates$es.R2conditional <- rsquareds$Conditional

        # format table nicely
        estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
        estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)

        outputList$results2 <- estimatesOutput

    }

    if (showEffectSizesTable) {

        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(estimates$term), effectSizes)

    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}

#' @export
summaryh.lme <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }

    # example output: b = −2.88, SE = 0.32, t(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(stats::coef(summary(model))) # get estimates and put in dataframe
    if (ncol(estimates) < 5 & (class(model)[1] == "merModLmerTest")) {
        return(message("lmerTest failed to compute p values; use summary() to check; try refitting model with lme() from nlme package"))
    }
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    colnames(estimates) <- c('estimate', 'std.error', 'df', 'statistic', 'p.value') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL

    # effect size r (Kashdan & Steger, 2006)
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

    # leave df as integers
    estimatesRound$df <- round(estimates$df)

    formattedOutput <- paste0("b = ", estimatesRound$estimate,
                              ", SE = ", estimatesRound$std.error,
                              ", t(", estimatesRound$df, ")",
                              " = ", estimatesRound$statistic,
                              ", p ", estimatesRound$p.value,
                              ", r = ", estimatesRound$es.r)

    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
    }

    formattedOutputDf <- data.table(term = as.character(estimates$term),
                                    results = as.character(formattedOutput))

    outputList <- list(results = formattedOutputDf)

    if (showTable) {

        # effect size semi partial R (Edwards et al., 2008)
        anovaModel <- data.frame(stats::anova(model))
        colnames(anovaModel) <- tolower(colnames(anovaModel))
        Fs <- anovaModel$f.value # F-values for each effect (marginal = type 3 SS with Satterthwaite (requires lmerTest package))
        numDF <- anovaModel$numdf #numerator DFs
        denDF <- anovaModel$dendf #denominator DFs
        semiPartialREffect <- (numDF / denDF * Fs) / (1 + (numDF / denDF * Fs)) # effect sizes
        semiPartialREffect <- data.frame(term = rownames(anovaModel), es.partR2 = semiPartialREffect, stringsAsFactors = F)
        for (i in 1:nrow(semiPartialREffect)) {
            estimates$es.partR2[grepl(semiPartialREffect$term[i], estimates$term)] <- semiPartialREffect[i, "es.partR2"]
        }

        # piecewiseSEM (Nakagawa & Schielzeth, 2013)
        # rsquareds <- sem.model.fits(model)
        # estimates$es.R2marginal <- rsquareds$Marginal
        # estimates$es.R2conditional <- rsquareds$Conditional

        # format table nicely
        estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
        estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)

        outputList$results2 <- estimatesOutput

    }

    if (showEffectSizesTable) {

        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(estimates$term), effectSizes)

    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}

#' @export
summaryh.lmerModLmerTest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }

    # example output: b = −2.88, SE = 0.32, t(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(stats::coef(summary(model))) # get estimates and put in dataframe
    if (ncol(estimates) < 5 & (class(model)[1] == "merModLmerTest")) {
        return(message("lmerTest failed to compute p values; use summary() to check; try refitting model with lme() from nlme package"))
    }
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    colnames(estimates) <- c('estimate', 'std.error', 'df', 'statistic', 'p.value') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL

    # effect size r (Kashdan & Steger, 2006)
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

    # leave df as integers
    estimatesRound$df <- round(estimates$df)

    formattedOutput <- paste0("b = ", estimatesRound$estimate,
                              ", SE = ", estimatesRound$std.error,
                              ", t(", estimatesRound$df, ")",
                              " = ", estimatesRound$statistic,
                              ", p ", estimatesRound$p.value,
                              ", r = ", estimatesRound$es.r)

    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
    }

    formattedOutputDf <- data.table(term = as.character(estimates$term),
                                    results = as.character(formattedOutput))

    outputList <- list(results = formattedOutputDf)

    if (showTable) {

        # effect size semi partial R (Edwards et al., 2008)
        anovaModel <- data.frame(stats::anova(model))
        colnames(anovaModel) <- tolower(colnames(anovaModel))
        Fs <- anovaModel$f.value # F-values for each effect (marginal = type 3 SS with Satterthwaite (requires lmerTest package))
        numDF <- anovaModel$numdf #numerator DFs
        denDF <- anovaModel$dendf #denominator DFs
        semiPartialREffect <- (numDF / denDF * Fs) / (1 + (numDF / denDF * Fs)) # effect sizes
        semiPartialREffect <- data.frame(term = rownames(anovaModel), es.partR2 = semiPartialREffect, stringsAsFactors = F)
        for (i in 1:nrow(semiPartialREffect)) {
            estimates$es.partR2[grepl(semiPartialREffect$term[i], estimates$term)] <- semiPartialREffect[i, "es.partR2"]
        }

        # piecewiseSEM (Nakagawa & Schielzeth, 2013)
        # rsquareds <- sem.model.fits(model)
        # estimates$es.R2marginal <- rsquareds$Marginal
        # estimates$es.R2conditional <- rsquareds$Conditional

        # format table nicely
        estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
        estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)

        outputList$results2 <- estimatesOutput

    }

    if (showEffectSizesTable) {

        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(estimates$term), effectSizes)

    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}

#' @export
summaryh.lmerTest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }

    # example output: b = −2.88, SE = 0.32, t(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(stats::coef(summary(model))) # get estimates and put in dataframe
    if (ncol(estimates) < 5 & (class(model)[1] == "merModLmerTest")) {
        return(message("lmerTest failed to compute p values; use summary() to check; try refitting model with lme() from nlme package"))
    }
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    colnames(estimates) <- c('estimate', 'std.error', 'df', 'statistic', 'p.value') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL

    # effect size r (Kashdan & Steger, 2006)
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

    # leave df as integers
    estimatesRound$df <- round(estimates$df)

    formattedOutput <- paste0("b = ", estimatesRound$estimate,
                              ", SE = ", estimatesRound$std.error,
                              ", t(", estimatesRound$df, ")",
                              " = ", estimatesRound$statistic,
                              ", p ", estimatesRound$p.value,
                              ", r = ", estimatesRound$es.r)

    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
    }

    formattedOutputDf <- data.table(term = as.character(estimates$term),
                                    results = as.character(formattedOutput))

    outputList <- list(results = formattedOutputDf)

    if (showTable) {

        # effect size semi partial R (Edwards et al., 2008)
        anovaModel <- data.frame(stats::anova(model))
        colnames(anovaModel) <- tolower(colnames(anovaModel))
        Fs <- anovaModel$f.value # F-values for each effect (marginal = type 3 SS with Satterthwaite (requires lmerTest package))
        numDF <- anovaModel$numdf #numerator DFs
        denDF <- anovaModel$dendf #denominator DFs
        semiPartialREffect <- (numDF / denDF * Fs) / (1 + (numDF / denDF * Fs)) # effect sizes
        semiPartialREffect <- data.frame(term = rownames(anovaModel), es.partR2 = semiPartialREffect, stringsAsFactors = F)
        for (i in 1:nrow(semiPartialREffect)) {
            estimates$es.partR2[grepl(semiPartialREffect$term[i], estimates$term)] <- semiPartialREffect[i, "es.partR2"]
        }

        # piecewiseSEM (Nakagawa & Schielzeth, 2013)
        # rsquareds <- sem.model.fits(model)
        # estimates$es.R2marginal <- rsquareds$Marginal
        # estimates$es.R2conditional <- rsquareds$Conditional

        # format table nicely
        estimatesOutput <- data.frame(lapply(estimates[, -1], round, decimal + 1))
        estimatesOutput <- data.table(term = as.character(estimates$term), estimatesOutput)

        outputList$results2 <- estimatesOutput

    }

    if (showEffectSizesTable) {

        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(estimates$term), effectSizes)

    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}

#' @export
summaryh.rma.uni <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, confInterval = NULL, ...) {

  # ensure significant digits with sprintf
  digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
  if (decimal <= 2) {
    pdigits <- paste0("%.", 3, "f")
  } else {
    pdigits <- paste0("%.", decimal, "f")
  }

  # example output: g = 0.12 (95% CI = [0.01, 1.01])
  estimates <- data.frame(estimate = as.numeric(model$b),
                          p.value = as.numeric(model$pval),
                          ciLower = as.numeric(model$ci.lb),
                          ciUpper = as.numeric(model$ci.ub)) # get estimates and put in dataframe

  # make a copy of estimates and convert to correct dp
  estimatesCopy <- estimates
  estimatesRound <- estimatesCopy
  estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
  estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
  estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
  estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
  # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

  # fix p values
  estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
  estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
  estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

  formattedOutput <- paste0("g = ", estimatesRound$estimate,
                            " (95% CI = [", estimatesRound$ciLower, ", ", estimatesRound$ciUpper, "])")

  if (estimates$p.value %between% c(0.01, 0.05)) {
    stars <- "*  "
  } else if (estimates$p.value %between% c(0.001, 0.01)) {
    stars <- "**"
  } else if (estimates$p.value < 0.001) {
    stars <- "***"
  } else {
    stars <- "   "
  }

  formattedOutput2 <- paste0(estimatesRound$estimate, ' ', stars, " [", estimatesRound$ciLower, ", ", estimatesRound$ciUpper, "]")

  # convert hyphens to minus (only possible on UNIX systems)
  if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
    formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
  }

  formattedOutputDf <- data.table(results = as.character(formattedOutput), results2 = as.character(formattedOutput2))

  outputList <- list(results = formattedOutputDf)

  if (showTable) {

    # format table nicely
    estimatesOutput <- data.frame(lapply(estimates, round, decimal + 1))
    estimatesOutput <- data.table(estimatesOutput)
    outputList$results2 <- estimatesOutput
  }

  if (showEffectSizesTable) {

    outputList$effectSizes <- NULL

  }

  options(scipen = 0) # enable scientific notation
  if (length(outputList) > 1) {
    return(outputList)
  } else {
    return(formattedOutputDf)
  }

}

# library(metaBMA)
# mf <- meta_fixed(towels$logOR, towels$SE, towels$study, d = "halfnorm",d.par = c(0, .3))
# mf
# summaryh(mf)
#' @export
summaryh.meta_fixed <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, confInterval = NULL, ...) {

  # ensure significant digits with sprintf
  digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
  if (decimal <= 2) {
    pdigits <- paste0("%.", 3, "f")
  } else {
    pdigits <- paste0("%.", decimal, "f")
  }

  # example output: d = 0.12 (95% HPD = [0.01, 1.01])
  estimatesTemp <- data.frame(model$estimates)
  estimates <- data.frame(estimate = as.numeric(estimatesTemp$Mean),
                          ciLower = as.numeric(estimatesTemp$HPD95lower),
                          ciUpper = as.numeric(estimatesTemp$HPD95upper),
                          bf = model$BF) # get estimates and put in dataframe

  # make a copy of estimates and convert to correct dp
  estimatesCopy <- estimates
  estimatesRound <- estimatesCopy
  estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
  estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
  estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
  estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
  # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

  # fix p values
  # estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
  # estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
  # estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
  #
  formattedOutput <- paste0("d = ", estimatesRound$estimate,
                            " (95% HPD = [", estimatesRound$ciLower, ", ", estimatesRound$ciUpper, "])")

  # if (estimates$p.value %between% c(0.01, 0.05)) {
  #     stars <- "*  "
  # } else if (estimates$p.value %between% c(0.001, 0.01)) {
  #     stars <- "**"
  # } else if (estimates$p.value < 0.001) {
  #     stars <- "***"
  # } else {
  #     stars <- "   "
  # }
  #

  formattedOutput2 <- paste0(estimatesRound$estimate, " [", estimatesRound$ciLower, ", ", estimatesRound$ciUpper, "]")

  # convert hyphens to minus (only possible on UNIX systems)
  if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
    formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
  }

  formattedOutputDf <- data.table(results = as.character(formattedOutput), results2 = as.character(formattedOutput2))

  outputList <- list(results = formattedOutputDf)

  if (showTable) {

    # format table nicely
    estimatesOutput <- data.frame(lapply(estimates, round, decimal + 1))
    estimatesOutput <- data.table(estimatesOutput)
    outputList$results2 <- estimatesOutput
  }

  if (showEffectSizesTable) {

    outputList$effectSizes <- NULL

  }

  options(scipen = 0) # enable scientific notation
  if (length(outputList) > 1) {
    return(outputList)
  } else {
    return(formattedOutputDf)
  }

}



# library(metaBMA)
# model <- meta_random(towels$logOR, towels$SE, towels$study, d = "halfnorm",d.par = c(0, .3))
# model
# class(model)
# summaryh(model)
#' @export
summaryh.meta_random <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, confInterval = NULL, ...) {

  # ensure significant digits with sprintf
  digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
  if (decimal <= 2) {
    pdigits <- paste0("%.", 3, "f")
  } else {
    pdigits <- paste0("%.", decimal, "f")
  }

  # example output: d = 0.12 (95% HPD = [0.01, 1.01])
  estimatesTemp <- data.frame(model$estimates)
  estimates <- data.frame(estimate = as.numeric(estimatesTemp$Mean),
                          ciLower = as.numeric(estimatesTemp$HPD95lower),
                          ciUpper = as.numeric(estimatesTemp$HPD95upper),
                          bf = model$BF) # get estimates and put in dataframe

  # make a copy of estimates and convert to correct dp
  estimatesCopy <- estimates
  estimatesRound <- estimatesCopy
  estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
  estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
  estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
  estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
  # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

  # fix p values
  # estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
  # estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
  # estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
  #
  formattedOutput <- paste0("d = ", estimatesRound$estimate,
                            " (95% HPD = [", estimatesRound$ciLower, ", ", estimatesRound$ciUpper, "])")

  # if (estimates$p.value %between% c(0.01, 0.05)) {
  #     stars <- "*  "
  # } else if (estimates$p.value %between% c(0.001, 0.01)) {
  #     stars <- "**"
  # } else if (estimates$p.value < 0.001) {
  #     stars <- "***"
  # } else {
  #     stars <- "   "
  # }
  #

  formattedOutput2 <- paste0(estimatesRound$estimate, " [", estimatesRound$ciLower, ", ", estimatesRound$ciUpper, "]")

  # convert hyphens to minus (only possible on UNIX systems)
  if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
    formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
  }

  formattedOutputDf <- data.table(results = as.character(formattedOutput), results2 = as.character(formattedOutput2))

  outputList <- list(results = formattedOutputDf)

  if (showTable) {

    # format table nicely
    estimatesOutput <- data.frame(lapply(estimates, round, decimal + 1))
    estimatesOutput <- data.table(estimatesOutput)
    outputList$results2 <- estimatesOutput
  }

  if (showEffectSizesTable) {

    outputList$effectSizes <- NULL

  }

  options(scipen = 0) # enable scientific notation
  if (length(outputList) > 1) {
    return(outputList)
  } else {
    return(formattedOutputDf)
  }

}

#' @export
summaryh.htest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE, ...) {
  if (grepl("t-test", model$method, ignore.case = T)) {
    reportTtest(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
  } else if (grepl("Pearson's product-moment correlation", model$method, ignore.case = T)) {
    reportCortestPearson(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
  } else if (grepl("Kendall's rank correlation tau", model$method, ignore.case = T)) {
    reportCortest(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
  } else if (grepl("Spearman's rank correlation rho", model$method, ignore.case = T)) {
    reportCortest(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
  } else if (grepl("Pearson's Chi-squared test", model$method, ignore.case = T)) {
    reportCHISQ(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
  }
}


#### specific htests
reportTtest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {

  # ensure significant digits with sprintf
  digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
  if (decimal <= 2) {
    pdigits <- paste0("%.", 3, "f")
  } else {
    pdigits <- paste0("%.", decimal, "f")
  }

  # example output: t(30) = 5.82, p < .001
  estimates <- data.frame(df = model$parameter,
                          statistic = model$statistic,
                          p.value = model$p.value)
  rownames(estimates) <- NULL

  # effect sizes
  estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
  estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d

  # make a copy of estimates and convert to correct dp
  estimatesCopy <- estimates[, -1]
  estimatesRound <- estimatesCopy
  estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
  estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
  estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
  estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
  # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

  # fix p values
  estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
  estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
  estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

  # leave df as integers
  estimatesRound$df <- round(estimates$df)

  formattedOutput <- paste0("t(", estimatesRound$df, ")",
                            " = ", estimatesRound$statistic,
                            ", p ", estimatesRound$p.value,
                            ", r = ", estimatesRound$es.r)

  # convert hyphens to minus (only possible on UNIX systems)
  if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
    formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
  }

  formattedOutputDf <- data.table(results = as.character(formattedOutput))

  outputList <- list(results = formattedOutputDf)

  if (showTable) {

    # format table nicely
    estimatesOutput <- data.frame(lapply(estimates, round, decimal + 1))
    estimatesOutput <- data.table(term = as.character(model$data.name), estimatesOutput)
    outputList$results2 <- estimatesOutput
  }

  if (showEffectSizesTable) {

    # get all other effect sizes
    effectSizes <- es(r = abs(round(as.numeric(estimates$es.r), decimal + 1)), decimal = decimal, msg = F)
    outputList$effectSizes <- data.table(term = as.character(model$data.name), effectSizes)

  }
  options(scipen = 0) # enable scientific notation
  if (length(outputList) > 1) {
    return(outputList)
  } else {
    return(formattedOutputDf)
  }

}

#' @importFrom compute.es chies
reportCHISQ <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {

  # ensure significant digits with sprintf
  digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
  if (decimal <= 2) {
    pdigits <- paste0("%.", 3, "f")
  } else {
    pdigits <- paste0("%.", decimal, "f")
  }

  # example output: r(30) = 0.82, p < .001
  estimates <- data.frame(df = model$parameter,
                          statistic = model$statistic,
                          p.value = model$p.value)
  rownames(estimates) <- NULL

  # effect sizes
  computeES <- chies(chi.sq = estimates$statistic, n = sum(model$observed), verbose = F)
  estimates$es.r <- computeES$r

  # make a copy of estimates and convert to correct dp
  estimatesCopy <- estimates[, -1]
  estimatesRound <- estimatesCopy
  estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
  estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
  estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
  estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
  # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

  # fix p values
  estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
  estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
  estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

  # leave df as integers
  estimatesRound$df <- round(estimates$df)

  formattedOutput <- paste0("X2(", estimatesRound$df, ")",
                            " = ", estimatesRound$statistic,
                            ", p ", estimatesRound$p.value,
                            ", r = ", estimatesRound$es.r)

  # convert hyphens to minus (only possible on UNIX systems)
  if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
    formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
  }

  formattedOutputDf <- data.table(results = as.character(formattedOutput))

  outputList <- list(results = formattedOutputDf)

  if (showTable) {

    # format table nicely
    estimatesOutput <- data.frame(lapply(estimates, round, decimal + 1))
    estimatesOutput <- data.table(term = as.character(model$data.name), estimatesOutput)
    outputList$results2 <- estimatesOutput
  }

  if (showEffectSizesTable) {

    # get all other effect sizes
    effectSizes <- es(r = abs(round(as.numeric(estimates$es.r), decimal + 1)), decimal = decimal, msg = F)
    outputList$effectSizes <- data.table(term = as.character(model$data.name), effectSizes)

  }
  options(scipen = 0) # enable scientific notation
  if (length(outputList) > 1) {
    return(outputList)
  } else {
    return(formattedOutputDf)
  }

}

reportCortest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {

  # ensure significant digits with sprintf
  digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
  if (decimal <= 2) {
    pdigits <- paste0("%.", 3, "f")
  } else {
    pdigits <- paste0("%.", decimal, "f")
  }

  # example output: r = 0.82, p < .001
  estimates <- data.frame(estimate = model$estimate,
                          p.value = model$p.value,
                          statistic = model$statistic)

  rownames(estimates) <- NULL

  # effect sizes
  estimates$es.d <- es(r = estimates$estimate, msg = F, decimal = decimal)$d # d

  # make a copy of estimates and convert to correct dp
  estimatesCopy <- estimates
  estimatesRound <- estimatesCopy
  estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
  estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
  estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
  estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
  # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

  # fix p values
  estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
  estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
  estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

  formattedOutput <- paste0("r = ", estimatesRound$estimate,
                            ", p ", estimatesRound$p.value)

  # convert hyphens to minus (only possible on UNIX systems)
  if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
    formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
  }

  formattedOutputDf <- data.table(results = as.character(formattedOutput))

  outputList <- list(results = formattedOutputDf)

  if (showTable) {

    # format table nicely
    # estimatesOutput <- apply(estimates, 2, round, decimal + 1)
    estimatesOutput <- data.frame(lapply(estimates, round, decimal + 1))
    estimatesOutput <- data.table(term = as.character(model$data.name), estimatesOutput)
    outputList$results2 <- estimatesOutput
  }

  if (showEffectSizesTable) {

    # get all other effect sizes
    effectSizes <- es(r = abs(round(as.numeric(estimates$estimate), decimal + 1)), decimal = decimal, msg = F)
    outputList$effectSizes <- data.table(term = as.character(model$data.name), effectSizes)

  }
  options(scipen = 0) # enable scientific notation
  if (length(outputList) > 1) {
    return(outputList)
  } else {
    return(formattedOutputDf)
  }

}

reportCortestPearson <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {

  # ensure significant digits with sprintf
  digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
  if (decimal <= 2) {
    pdigits <- paste0("%.", 3, "f")
  } else {
    pdigits <- paste0("%.", decimal, "f")
  }

  # example output: r(30) = 0.82, p < .001

  estimates <- data.frame(df = model$parameter,
                          estimate = model$estimate,
                          estimateLow = model$conf.int[1],
                          estimateUpper = model$conf.int[2],
                          p.value = model$p.value,
                          statistic = model$statistic)

  rownames(estimates) <- NULL

  # effect sizes
  estimates$es.d <- (2 * estimates$statistic) / sqrt(estimates$df) # d

  # make a copy of estimates and convert to correct dp
  estimatesCopy <- estimates[, -1]
  estimatesRound <- estimatesCopy
  estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
  estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
  estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
  estimatesRound[abs(estimatesCopy) < 0.0000001] <- 0
  # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])

  # fix p values
  estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
  estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
  estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)

  # leave df as integers
  estimatesRound$df <- round(estimates$df)

  formattedOutput <- paste0("r(", estimatesRound$df, ")",
                            " = ", estimatesRound$estimate,
                            ", p ", estimatesRound$p.value)

  # convert hyphens to minus (only possible on UNIX systems)
  if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
    formattedOutput <- gsub("-", replacement = "\u2212", formattedOutput)
  }

  formattedOutputDf <- data.table(results = as.character(formattedOutput))

  outputList <- list(results = formattedOutputDf)

  if (showTable) {

    # format table nicely
    estimatesOutput <- data.frame(lapply(estimates, round, decimal + 1))
    estimatesOutput <- data.table(term = as.character(model$data.name), estimatesOutput)
    outputList$results2 <- estimatesOutput
  }

  if (showEffectSizesTable) {

    # get all other effect sizes
    effectSizes <- es(r = abs(round(as.numeric(estimates$estimate), decimal + 1)), decimal = decimal, msg = F)
    outputList$effectSizes <- data.table(term = as.character(model$data.name), effectSizes)

  }
  options(scipen = 0) # enable scientific notation
  if (length(outputList) > 1) {
    return(outputList)
  } else {
    return(formattedOutputDf)
  }

}
