#' @title Identify outliers using robust median absolute deviation approach
#' @name outliersMAD
#'
#' @description outliersMAD is used to identify outliers in vectors using Leys et al.'s (2003) median absolute deviation approach.
#'
#' @param x a vector of numbers
#' @param MADCutOff value to use as cutoff (Leys e tal. recommend 2.5 or 3.0 as default)
#' @param replaceOutliersWith if value is an outlier, what to replace it with? NA by default
#' @param showMADValues if TRUE, will show deviation score of each value
#' @param outlierIndices return index/position of outlier
#' @param bConstant a constant linked to the assumption of normality of the data, disregarding the abnormality induced by outliers
#' @param digits how many digits/decimals to round output to
#'
#' @return A vector with outliers identified (default converts outliers to NA)
#'
#' @details We can identify and remove outliers in our data by identifying data points that are too extreme—either too many standard deviations (SD) away from the mean or too many median absolute deviations (MAD) away from the median. The SD approach might not be ideal with extreme outliers, whereas the MAD approach is much more robust (for comparison of both approaches, see Leys et al., 2013, Journal of Experimental Social Psychology).
#'
#' @references \itemize{
#' \item Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. Journal of Experimental Social Psychology, 49(4), 764-766. doi:10.1016/j.jesp.2013.03.013 (\url{https://www.sciencedirect.com/science/article/pii/S0022103113000668})}
#' @seealso \code{\link{outliersZ}}
#'
#' @author Hause Lin
#'
#' @export
#'
#' @usage
#' outliersMAD(x, MADCutOff = 2.5, replaceOutliersWith = NA,
#' showMADValues = FALSE, outlierIndices = FALSE, bConstant = 1.4826, digits = 2)
#'
#' @examples
#' example <- c(1, 3, 3, 6, 8, 10, 10, 1000) # 1000 is an outlier
#' outliersMAD(example)
#' outliersMAD(example, MADCutOff = 3.0)
#' outliersMAD(example, MADCutOff = 2.5, replaceOutliersWith = -999)
#' outliersMAD(example, MADCutOff = 1.5, outlierIndices = TRUE)
#' outliersMAD(example, MADCutOff = 1.5, showMADValues = TRUE)
outliersMAD <- function(x, MADCutOff = 2.5, replaceOutliersWith = NA, showMADValues = FALSE, outlierIndices = FALSE, bConstant = 1.4826, digits = 2) {
  # bConstant: usually, b = 1.4826, a constant linked to the assumption of normality of the data, disregarding the abnormality induced by out- liers (Rousseeuw & Croux, 1993).

  # compute number of absolute MADs away for each value: formula: abs( ( x - median(x) ) )/ mad(x)
  absMADAway <- abs((x - stats::median(x, na.rm = T))/stats::mad(x, constant = bConstant, na.rm = T))
  # subset data that has absMADAway greater than the MADCutOff and replace them with replace
  x[absMADAway > MADCutOff] <- replaceOutliersWith
  outliers <- length(x[absMADAway > MADCutOff])
  if (showMADValues) { # if values == TRUE, return number of mads for each value
    message("Showing absolute MAD from median for each value.")
    message(paste0(outliers, " outliers detected."))
    return(round(absMADAway, digits))
  } else if (outlierIndices) {
    message("Showing indices of outliers.")
    return(which(is.na(x)))
  } else {
    message(paste0(outliers, " outliers detected."))
    message(paste0("Outliers replaced with ", replaceOutliersWith))
    return(round(x, digits)) # otherwise, return original with outliers replaced
  }
}
