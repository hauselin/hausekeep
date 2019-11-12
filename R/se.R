#' @title Compute standard errors (between-subjects)
#' @name stderror
#' @description se is used to compute the standard error(s) for one or more variables, for one or more groups in a dataframe.
#' @param data a dataframe
#' @param measurevar the name(s) of column(s) that contain the variable to be summariezed
#' @param groupvars a vector containing names of columns that contain grouping variables
#' @param na.rm boolean that indicates whether to ignore NA values
#' @param conf.interval confidence interval range
#' @param tonumeric whether to convert variables/columns to numeric whenever possible
#'
#' @usage
#' stderror(data = NULL, measurevar, groupvars = NULL,
#' na.rm = TRUE, conf.interval = 0.95, tonumeric = TRUE)
#' @note Code adapted from R cookbook (see references)
#' @author Hause Lin
#' @import data.table
#' @importFrom dplyr tbl_df
#' @importFrom stats sd
#' @return A dataframe (tibble)
#' @examples
#' stderror(data = mtcars, measurevar = c("mpg", "disp"), groupvars = c("cyl", "am"))
#' stderror(data = ChickWeight, measurevar = "weight", groupvars = "Diet")
#' @export
stderror <- function(data = NULL, measurevar, groupvars = NULL, na.rm = TRUE, conf.interval = 0.95, tonumeric = TRUE) {

  # convert to datatable and tibble
  data <- data.table::data.table(data)

  # function to compute N without NAs
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm)
      sum(!is.na(x))
    else length(x)
  }

  resultsList <- list() # empty list to store results

  for (i in 1:length(measurevar)) {

    if (sum( data.frame(data)[, measurevar[i]] %in% c(Inf, -Inf)) > 0) { # if measurvar contains Inf or -Inf, stop the script
      stop(paste0("\nInf or -Inf is in ", measurevar[i], " variable"))
    }

    # compute mean by group
    datac <- data[, .(unlist(lapply(.SD, length2, na.rm = na.rm)),
                      unlist(lapply(.SD, mean, na.rm = na.rm)),
                      unlist(lapply(.SD, sd, na.rm = na.rm))),
                  by = groupvars, .SDcols = measurevar[i]]

    setnames(datac, c(groupvars, "N", measurevar[i], "sd")) # rename column names

    setkeyv(datac, groupvars) # sort table

    datac[, se := sd / sqrt(N)] # compute standard error

    ciMult <- stats::qt(conf.interval / 2 + 0.5, unlist(datac$N) - 1)
    datac[, ci := se * ciMult]

    if (tonumeric) {
      # convert columns to numeric class if possible, else, leave as character
      oldwarning <- getOption("warn")
      options(warn = -1)
      for (j in 1:(ncol(datac)-4)) { # exclude last few columns (outcome, sd, se, ci)
        if (sum(is.na(as.numeric(as.character(datac[[j]])))) == 0) {
          datac[[j]] <- as.numeric(datac[[j]])
        } else {
          datac[[j]] <- as.character(datac[[j]])
        }
      }
      options(warn = oldwarning)
    }

    resultsList[[measurevar[i]]] <- tbl_df(datac)

  }

  if (length(measurevar) == 1) {
    return(resultsList[[measurevar[1]]])
  } else {
    return(resultsList)
  }

}


#' @title Norms data within groups
#' @name normWithin
#' @param data a dataframe
#' @param idvar the name of the column that identifies each subject
#' @param measurevar the name of a column that contains the variable to be summariezed
#' @param betweenvars a vector containing names of columns that are between-subjects variables
#' @param na.rm boolean that indicates whether to ignore NA's
#' @return a data.frame
#' @description normWithin norms the data within specified groups in a data frame
#' @importFrom dplyr left_join
#' @import data.table
#' @examples
#' \dontrun{
#' normWithin(data = sleep, idvar = "ID", measurevar = "extra", betweenvars = "group")
#' }
normWithin <- function (data = NULL, idvar, measurevar, betweenvars = NULL, na.rm = TRUE) {
  # Norms the data within specified groups in a data frame; it normalizes each
  # subject (identified by idvar) so that they have the same mean, within each group
  # specified by betweenvars.
  # norm data (this function will only be used by seWithin, and won't have to be called directly)

  data <- data.table::data.table(data)
  setkeyv(data, idvar) # sort by idvar

  data.subjMean <- data[, .(unlist(lapply(.SD, mean, na.rm = na.rm))), by = c(idvar, betweenvars), .SDcols = measurevar] # compute mean for each subject
  setnames(data.subjMean, c(idvar, betweenvars,'subjMean'))
  dataNew <- left_join(data, data.subjMean)
  dataNew <- data.table::data.table(dataNew)
  setkeyv(dataNew, c(idvar, betweenvars)) # sort

  measureNormedVar <- paste0(measurevar, "Normed")

  # dataNew <- data.frame(dataNew)
  # dataNew[, measureNormedVar] <- dataNew[, measurevar] - unlist(data[, "subjMean"]) + mean(data[, measurevar], na.rm = na.rm)

  dataNew[, (measureNormedVar) := get(measurevar) - subjMean + mean(get(measurevar), na.rm = T)]

  dataNew$subjMean <- NULL

  return(data.frame(dataNew))
}



#' @title Compute standard errors (within-subjects)
#' @name seWithin
#'
#' @description seWithin is used to compute the standard error(s) for one or more variables, for one or more groups in a dataframe, handling within-subjects variables by removing inter-subject variability
#'
#' @param data a dataframe
#' @param measurevar the name(s) of column(s) that contain the variable to be summariezed
#' @param betweenvars a vector containing names of columns that are between-subjects variables
#' @param withinvars a vector containing names of columns that are within-subjects variables
#' @param idvar the name of the column that identifies each subject
#' @param na.rm boolean that indicates whether to ignore NA's
#' @param conf.interval confidence interval range
#' @param shownormed whether to show noramlized results
#'
#' @seealso \code{\link{stderror}}
#'
#' @author Hause Lin
#'
#' @note
#' Code adapted from R cookbook (see references).
#'
#' @references
#' \url{http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/}
#'
#' @return A data.table
#'
#' @import data.table
#' @importFrom dplyr tbl_df
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate_if
#' @importFrom stats sd
#' @export
#' @seealso \code{\link{stderror}}
#' @usage
#' seWithin(data = NULL, measurevar, betweenvars = NULL, withinvars = NULL,
#' idvar = NULL, na.rm = TRUE, conf.interval = 0.95, shownormed = FALSE)
#'
#' @examples
#' result <- seWithin(data = ChickWeight, measurevar = "weight",
#' betweenvars = "Diet", withinvars = "Time", idvar = "Chick")
#'
#' # multiple outcome variables
#' ChickWeight2 <- ChickWeight
#' ChickWeight2$weight2 <- ChickWeight2$weight * 100 # create a new outcome variable
#' result <- seWithin(data = ChickWeight2, measurevar = c("weight", "weight2"),
#' betweenvars = "Diet", withinvars = "Time", idvar = "Chick")
seWithin <- function (data = NULL, measurevar, betweenvars = NULL, withinvars = NULL, idvar = NULL, na.rm = TRUE, conf.interval = 0.95, shownormed = FALSE) {
  # within-subjects CI (normed and un-normed versions)
  ## Summarizes data, handling within-subjects variables by removing inter-subject variability.
  ## It will still work if there are no within-S variables.
  ## Gives count, un-normed mean, normed mean (with same between-group mean),
  ##   standard deviation, standard error of the mean, and confidence interval.
  ## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
  ##   data: a data frame.
  ##   measurevar: the name of a column that contains the variable to be summariezed
  ##   betweenvars: a vector containing names of columns that are between-subjects variables
  ##   withinvars: a vector containing names of columns that are within-subjects variables
  ##   idvar: the name of a column that identifies each subject (or matched subjects)
  ##   na.rm: a boolean that indicates whether to ignore NA's
  ##   conf.interval: the percent range of the confidence interval (default is 95%)
  ##    shownormed: whether to show the normed version of the outcome variable

  data <- data.frame(data) # convert to data.frame

  # Check if betweenvars and withinvars are factors
  factorvars <- sapply(data[, c(betweenvars, withinvars), drop = FALSE], FUN = is.factor)
  # Ensure that the betweenvars and withinvars are factors
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ", paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }

  resultsList <- list() # empty list to store results

  for (i in 1:length(measurevar)) {

    # if measurvar contains Inf or -Inf, stop the script
    if (sum( data[, measurevar[i]] %in% c(Inf, -Inf)) > 0) {
      stop(paste0("\nInf or -Inf is in ", measurevar[i], " variable"))
    }

    # Get the means from the un-normed data
    datac <- stderror(data, measurevar[i], groupvars = c(betweenvars, withinvars), na.rm = na.rm, conf.interval = conf.interval, tonumeric = FALSE)

    # Drop all the unused columns (these will be calculated with normed data)
    datac$sd <- NULL
    datac$se <- NULL
    datac$ci <- NULL

    # Norm each subject's data
    ndata <- normWithin(data, idvar, measurevar[i], betweenvars, na.rm)

    # This is the name of the new column
    measurevar_n <- paste(measurevar[i], "Normed", sep = "")

    # Collapse the normed data - now we can treat between and within vars the same
    ndatac <- stderror(ndata, measurevar_n, groupvars = c(betweenvars, withinvars), na.rm = na.rm, conf.interval = conf.interval, tonumeric = FALSE)

    # Apply correction from Morey (2008) to the standard error and confidence interval
    # Get the product of the number of conditions of within-S variables
    nWithinGroups <- prod(vapply(ndatac[,withinvars, drop = FALSE], FUN = function(x) length(levels(x)), FUN.VALUE = numeric(1)))
    correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )

    ndatacTbl <- data.table::data.table(ndatac)

    # Apply the correction factor
    # setnames(ndatacTbl, c("sd", "se"), c("stdev", "stderror"))
    # print(ndatacTbl)
    ndatacTbl[, `:=` (sd = sd * correctionFactor, se = se * correctionFactor, ci = ci * correctionFactor)]
    # print(ndatacTbl)
    # setnames(ndatacTbl, c("stdev", "stderror"), c("sd", "se"))

    # Combine the un-normed means with the normed results
    merged <- left_join(datac, ndatacTbl)
    merged <- mutate_if(merged, is.factor, as.character) #if factor, convert to character
    merged[order( unlist((merged[, 1])), decreasing =  F), ] #arrange by first column
    merged <- data.table::data.table(merged)
    message("Factors have been converted to characters.")

    # convert columns to numeric class if possible, else, leave as character
    oldwarning <- getOption("warn")
    options(warn = -1)
    for (j in 1:(ncol(merged)-4)) { # exclude last few columns (outcome, sd, se, ci)
      if (sum(is.na(as.numeric(as.character(merged[[j]])))) == 0) {
        merged[[j]] <- as.numeric(merged[[j]])
      } else {
        merged[[j]] <- as.character(merged[[j]])
      }
    }
    options(warn = oldwarning)

    # whether to show normed version
    if (shownormed == FALSE) {
      # print(measurevar_n)
      merged[, (measurevar_n) := NULL]
    }

    resultsList[[measurevar[i]]] <- data.table(merged)
    message(cat("Confidence intervals: ", conf.interval, sep = ""))

  }

  if (length(measurevar) == 1) {
    print(resultsList[[measurevar[1]]])
    return(resultsList[[measurevar[1]]])
  } else {
    print(resultsList)
    return(resultsList)
  }

}
