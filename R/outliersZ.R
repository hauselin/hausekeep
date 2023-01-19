#' @title Identify outliers using Z-score cut-off method.
#' @name outliersZ
#'
#' @description outliersZ is used to identify outliers in vectors using Z-score cut-off
#'
#' @param x a vector of numbers
#' @param zCutOff value to use as cutoff (1.96 is a common value)
#' @param replaceOutliersWith if value is an outlier, what to replace it with? NA by default
#' @param outlierIndices return index/position of outlier
#' @param showZValues if TRUE, will show z score of each value
#' @param digits how many digits/decimals to round output to
#'
#' @return A vector with outliers identified (default converts outliers to NA)
#'
#' @note This detection method is not as robust as the median absolute deviation outlier detection method.
#'
#' @seealso \code{\link{outliers_mad}}
#'
#' @author Hause Lin
#'
#' @export
#'
#' @usage
#' outliersZ(x, zCutOff = 1.96, replaceOutliersWith = NA,
#' outlierIndices = FALSE, showZValues = FALSE, digits = 2)
#'
#' @examples
#' example <- c(1, 3, 3, 6, 8, 10, 10, 1000) # 1000 is an outlier
#' outliersZ(example)
#' outliersZ(example, zCutOff = 3.0)
#' outliersZ(example, zCutOff = 1.0, replaceOutliersWith = -999)
#' outliersZ(example, zCutOff = 1.0, outlierIndices = TRUE)
#' outliersZ(example, zCutOff = 1.0, showZValues = TRUE)
outliersZ <- function(x, zCutOff = 1.96, replaceOutliersWith = NA, outlierIndices = FALSE, showZValues = FALSE, digits = 2) {

  # compute standard deviation (sample version n = n [not n-1])
  stdev <- sqrt(sum((x - mean(x, na.rm = T))^2, na.rm = T) / sum(!is.na(x)))
  # compute Z values for each value
  Zvals <- (x - mean(x, na.rm = T)) / stdev
  absZ <- abs(Zvals)
  # subset data that has absZ greater than the zCutOff and replace them with replace
  # can also replace with other values (such as max/mean of data)
  x[absZ > zCutOff] <- replaceOutliersWith
  outliers <- length(x[absZ > zCutOff])

  if (showZValues) {
    # message("Showing absolute z-scores for each value.")
    # message(paste0(outliers, " outliers detected."))
    return(round(Zvals, digits)) # if values == TRUE, return z score for each value
  } else if (outlierIndices) {
    message("Showing indices of outliers.")
    return(which(is.na(x)))
  } else {
    # message(paste0(outliers, " outliers detected."))
    # message(paste0("Outliers replaced with ", replaceOutliersWith))
    return(round(x, digits)) # otherwise, return values with outliers replaced
  }
}
