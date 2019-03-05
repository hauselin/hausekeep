#' @title Fit Wagenmaker et al.'s (2007) EZ-diffusion model to multiple subjects and groups using raw single-trial data.
#' @name fit_ezddm
#'
#' @description Fits Wagenmaker et al.'s (2007) EZ-diffusion model for two-choice response time tasks. To use the function, ensure your dataframe is in long form, has single-trial reaction time (in seconds) and responses (coded as 0 or 1) on each row. You can use the function to fit the EZ-diffusion model to just a single subject or multiple subjects, and separately for each experimental condition (see below for examples).
#'
#' @param data data object with reaction time and accuracy variables (long form data expected)
#' @param rts specify in characters the name of the reaction time column (reaction time must be in seconds)
#' @param responses specify in characters the name of the accuracy column (coded as 0/1)
#' @param id specify in characters the name of your subject/id column (if not specified, assumes data [all rows] belong to a single subject)
#' @param group specify in characters the name of your column(s) indicating various conditions (default = NULL)
#' @param simCheck simulate data (n = 1000) with estimated parameters (using rdiffusion from rtdists package) to check model fit (default = FALSE)
#' @param decimal round parameter estimates (default = 4)
#'
#' @return A datatable
#'
#' @details Details to follow
#'
#' @note Notes to follow
#'
#' @references \itemize{
#' \item Wagenmakers, E. J., van der Maas, H. L., & Grasman, R. P. (2007). An EZ-diffusion model for response time and accuracy. Psychonomic Bulletin & Review, 14(1), 3-22. doi:10.3758/BF03194023 (\url{https://link.springer.com/article/10.3758/BF03194023})
#' }
#'
#' @author Hause Lin
#'
#' @seealso \code{\link{ezddm}}
#'
#' @import rtdists
#' @import data.table
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom dtplyr tbl_dt
#' @importFrom dplyr %>%
#' @export
#'
#' @usage
#' fit_ezddm(data, rts, responses, id = NULL, group = NULL, simCheck = FALSE, decimal = 4)
#'
#' @examples
#' library(rtdists) # load package to simulate data with diffusion parameters
#' data1 <- rdiffusion(n = 100, a = 2, v = 0.3, t0 = 0.5, z = 0.5 * 2) # simulate data
#' data2 <- rdiffusion(n = 100, a = 2, v = -0.3, t0 = 0.5, z = 0.5 * 2) # simulate data
#' dataAll <- rbind(data1, data2) # join data
#'
#' dataAll$response <- ifelse(dataAll$response == "upper", 1, 0) # convert responses to 1 and 0
#' dataAll$subject <- rep(c(1, 2), each = 100) # assign subject id
#'
#' dataAll$cond1 <- base::sample(c("a", "b"), 200, replace = TRUE) # randomly assign conditions a/b
#' dataAll$cond2 <- base::sample(c("y", "z"), 200, replace = TRUE) # randomly assign conditions y/z
#'
#' # fit model to just entire data set (assumes all data came from 1 subject)
#' fit_ezddm(data = dataAll, rts = "rt", responses = "response")
#' # fit model to just entire data set (assumes all data came from 1 subject)
#' fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject")
#' # fit model to each subject by cond1
#' fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = "cond1")
#' # fit model to each subject by cond1,cond2
#' fit_ezddm(data = dataAll, rts = "rt", responses = "response",
#' id = "subject", group = c("cond1", "cond2"))
fit_ezddm <- function(data, rts, responses, id = NULL, group = NULL, simCheck = FALSE, decimal = 4) {

  # message("Fits EZ-diffusion model (Wagenmaker et al., 2007, Psychonomic Bulletin & Review).\nResponses or choice must be coded as 0 (lower bound) or 1 (upper bound).")

  data <- tbl_dt(data)

  # create new variables
  data$rtCol <- data[, get(rts)]
  data$responseCol <- data[, get(responses)]

  # check if rt is in seconds
  if (mean(data$rtCol, na.rm = T) > 100) {
    stop("Check if reaction time is in seconds, not milliseconds!")
  }

  # if no id variable provided, assume it's just one subject's data
  if (is.null(id)) {
    id <- "temporary_subject"
    data$temporary_subject <- 1
  }

  # remove rts or responses rows
  data <- data[!is.na(rtCol), ]
  data <- data[!is.na(responseCol), ]

  # recode response accordingly (character and integer)
  if (data[, unique(responseCol)][1] %in% c(0, 1)) {
    data[, response_num := responseCol]
    data[, response_char := ifelse(response_num == 1, 'upper', 'lower')]
  } else if (data[, unique(responseCol)][1] %in% c('upper', 'lower')) {
    data[, response_char := as.character(responseCol)]
    data[, response_num := ifelse(response_char == "upper", 1, 0)]
  }

  # compute rt and responses
  behavOverall <- data[, .(response = round(mean(response_num, na.rm = T), 3),
                           rtOverall = round(mean(rtCol, na.rm = T), 3)),
                       by = c(id, group)]
  behav0 <- data[response_num == 0, .(rt0 = round(mean(rtCol, na.rm = T), 3)), by = c(id, group)]
  behav1 <- data[response_num == 1, .(rt1 = round(mean(rtCol, na.rm = T), 3)), by = c(id, group)]
  behav <- left_join(behavOverall, behav0, by = c(id, group))
  behav <- left_join(behav, behav1, by = c(id, group))

  # get grouping variables
  dataGroup <- data[, .(n = .N), by = c(id, group)]
  dataGroup0 <- data[response_num == 0, .(n0 = .N), by = c(id, group)]
  dataGroup1 <- data[response_num == 1, .(n1 = .N), by = c(id, group)]
  dataGroup <- left_join(dataGroup, dataGroup0, by = c(id, group))
  dataGroup <- left_join(dataGroup, dataGroup1, by = c(id, group))

  # for accurate responses (coded as 1), calculate mean RT and RT variance for each subject, each condition
  ddmRt <- data[response_num == 1, .(rt = mean(rtCol, na.rm = T), rtVar = stats::var(rtCol, na.rm = T)), by = c(id, group)]

  # calculate responses for each subject, each condition
  ddmAcc <- tbl_dt(data[, .(acc = mean(response_num, na.rm = T), n = .N), by = c(id, group)])

  if (sum(ddmAcc[, acc] %in% c(0.5, 1)) > 0) {
    n_corrected <- sum(ddmAcc[, acc] %in% c(0.5, 1))
    message(paste0("Mean accuracies (n = ", n_corrected, ") that are 0.5, or 1 have been adjusted slightly for model fitting."))
    ddmAcc[, accAdjust := 0]
    ddmAcc[acc %in% c(0.5, 1), accAdjust := 1]
  }

  # if acc is 1, apply edge correction
  ddmAcc[acc == 1, acc := edgeCorrect(n)] # edge correction
  # if acc is 0 or 50, add 0.001 to acc a bit so model fitting works
  ddmAcc[acc %in% c(0.5), acc := acc + 0.00001]

  dataForDDM <- left_join(ddmRt, ddmAcc, by = c(id, group))
  setDT(dataForDDM)

  # fit ez ddm model to each subject, each condition
  ddmResults <- dataForDDM[, ezddm(propCorrect = acc, rtCorrectVariance_seconds = rtVar, rtCorrectMean_seconds = rt),
                           by = c(id, group)]
  resultsFinal <- left_join(dataGroup, ddmResults, by = c(id, group))
  resultsFinal <- left_join(resultsFinal, ddmRt, by = c(id, group))
  resultsFinal <- left_join(resultsFinal, ddmAcc, by = c(id, group, "n"))
  resultsFinal <- left_join(resultsFinal, behav, by = c(id, group))
  resultsFinal <- dplyr::select(resultsFinal, -rt, -acc)

  setDT(resultsFinal) # ensure it's data table format
  setnames(resultsFinal, c("Ter", "rtVar"), c("t0_Ter", "rt1Var"))

  # simulate data
  if (simCheck) {
    resultsToSimulate <- copy(resultsFinal)
    simulatedData <- resultsToSimulate[a > 0 & t0_Ter > 0, rdiffusion(n = 1000, a = a * 10, v = v * 10, t0 = t0_Ter, z = 0.5 * a * 10), by = c(id, group)]
    simulatedData[, response_num := as.numeric()]
    simulatedData[response == 'upper', response_num := 1]
    simulatedData[response == 'lower', response_num := 0]
    simulateBehavOverall <- simulatedData[, .(responseSim = round(mean(response_num, na.rm = T), 3),
                                              rtOverallSim = round(mean(rt, na.rm = T), 3)),
                                          by = c(id, group)]
    simulateBehav0 <- simulatedData[response_num == 0, .(rt0Sim = round(mean(rt, na.rm = T), 3)), by = c(id, group)]
    simulateBehav1 <- simulatedData[response_num == 1, .(rt1Sim = round(mean(rt, na.rm = T), 3)), by = c(id, group)]
    simulateBehav <- left_join(simulateBehavOverall, simulateBehav0, by = c(id, group))
    simulateBehav <- left_join(simulateBehav, simulateBehav1, by = c(id, group))
    resultsFinal <- left_join(resultsFinal, simulateBehav, by = c(id, group))
    resultsFinal <- dplyr::select(resultsFinal, 1, group, n:rt1Var, response, responseSim, rtOverall, rtOverallSim, rt0, rt0Sim, rt1, rt1Sim)
  }

  # round results
  resultsFinal[, a := round(a, decimal)]
  resultsFinal[, v := round(v, decimal)]
  resultsFinal[, t0_Ter := round(t0_Ter, decimal)]
  resultsFinal[, rt1Var := round(rt1Var, 3)]

  # remove temporary_subject variable
  if (id == 'temporary_subject') {
    resultsFinal$temporary_subject <- NULL
  }

  return(tbl_dt(resultsFinal))

}


#' @title Fit Wagenmaker et al.'s (2007) EZ-diffusion model to single subject's aggregate/averaged data.
#' @name ezddm
#'
#' @description ezddm fits Wagenmaker et al.'s (2007) EZ-diffusion model for two-choice response time tasks using proportion correct, variance of correct response reaction times, and mean of correct reaction times (in seconds).
#'
#' @param propCorrect proportion correct (apply edge correction if necessary)
#' @param rtCorrectVariance_seconds variance of correct reaction times (in seconds)
#' @param rtCorrectMean_seconds mean of correct reaction times (in seconds)
#' @param nTrials number of trials (useful for edge correction) (optional)
#'
#' @return A dataframe with a (boundary parameter), v (drift rate parameter), and Ter (non-decision time parameter)
#'
#' @author Hause Lin
#'
#' @references \itemize{
#' \item Wagenmakers, E. J., van der Maas, H. L., & Grasman, R. P. (2007). An EZ-diffusion model for response time and accuracy. Psychonomic Bulletin & Review, 14(1), 3-22. doi:10.3758/BF03194023 (\url{https://link.springer.com/article/10.3758/BF03194023})
#' }
#'
#' @export
#' @examples
#' ezddm(.802, .112, .723)
ezddm <- function(propCorrect, rtCorrectVariance_seconds, rtCorrectMean_seconds, nTrials = NULL) {

    s <- 0.1 # s is scaling parameter (defaults to 0.1 in Ratcliff's models)
    s2 <- s^2 # variance

    v <- as.numeric(NA)
    a <- as.numeric(NA)
    Ter <- as.numeric(NA)

    # if propCorrect equals 0, 0.5, or 1, this method will not work, and an edge correction is required
    if (propCorrect %in% c(0, 0.5, 1)) {

        if (propCorrect == 0) {
            return(cat("Oops, propCorrect == 0. Can't fit model! D:"))
        } else if (propCorrect == 0.5) {
            cat("Oops, propCorrect == 0.5 (chance performance; drift will be close to 0). Added 0.00001 to propCorrect.\n")
            propCorrect <- propCorrect + 0.00001
        } else if (propCorrect == 1) {
            if (!is.null(nTrials)) {
                cat("Oops, propCorrect == 1. Applied edge correction.\n")
                propCorrect <- 1 - (1 / (2 * nTrials))
            } else {
                cat("Oops, propCorrect == 1. Edge correction required. Provide number of trials (nTrials).\n")
            }
        }

    }

    if (propCorrect != 1) {

        L <- stats::qlogis(propCorrect) # calculates logit
        x <- L * (L * propCorrect^2 - L * propCorrect + propCorrect - 0.5) / rtCorrectVariance_seconds
        v <- sign(propCorrect - 0.5) * s * x^(1/4) # drift rate
        a <- s2 * stats::qlogis(propCorrect)/v # threshold
        y <- -v*a/s2
        MDT <- (a/(2*v)) * (1-exp(y))/(1 + exp(y))
        Ter <- rtCorrectMean_seconds - MDT # non-decision time

    }

    return(data.frame(a, v, Ter))
}





edgeCorrect <- function(n) {
    return(1 - (1 / (2 * n))) # n: number of observations
}


# ezddm(.802, .112, .723)
# ezddm(.5, .112, .723)
# ezddm(.51, .112, .723)
# ezddm(0, .112, .723)
# ezddm(0.0001, .112, .723)
# ezddm(0, .112, .723)
# ezddm(0.005, .112, .723)
# ezddm(0.005, .112, .723)
# ezddm(1, .112, .723, 100)
#
# ezddm(0.8881988, 0.1005484, 0.9010186)
# library(EZ2)
# Data2EZ(.802, .112, .723)
# Data2EZ(.5, .112, .723)
# Data2EZ(0.8881988, 0.1005484, 0.9010186)
# Data2EZ(0.1, 0.1005484, 0.9010186)
# Data2EZ(0.00001, 0.1005484, 0.9010186)
# ezddm(0.000001, 0.1005484, 0.9010186)
# ezddm(0.00001, 0.1005484, 0.9010186)
# ezddm(0.5, 0.1005484, 0.9010186)
# ezddm(0.51, 0.1005484, 0.9010186)
# data.frame(Data2EZ(.802, .112, .723))
