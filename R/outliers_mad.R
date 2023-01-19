#' @title Identify outliers using robust median absolute deviation approach
#' @name outliers_mad
#'
#' @description outliers_mad is used to identify outliers in vectors using Leys et al.'s (2003) median absolute deviation approach.
#'
#' @param x a vector of numbers
#' @param threshold value to use as cutoff (Leys et al. recommend 2.5 or 3.0 as default)
#' @param replace_outlier_value if value is an outlier, what to replace it with? NA by default
#' @param show_mad_values if TRUE, will show deviation score of each value
#' @param show_outlier_indices if TRUE, return index/position of outliers
#' @param b_constant a constant linked to the assumption of normality of the data, disregarding the abnormality induced by outliers
#' @param digits how many digits to round output to
#' @param debug if TRUE, print messages (FALSE by default)
#'
#' @return A vector with outliers identified (default converts outliers to NA)
#'
#' @details We can identify and remove outliers in our data by identifying data points that are too extreme—either too many standard deviations (SD) away from the mean or too many median absolute deviations (MAD) away from the median. The SD approach might not be ideal with extreme outliers, whereas the MAD approach is much more robust (for comparison of both approaches, see Leys et al., 2013, Journal of Experimental Social Psychology).
#'
#' b_constant is usually 1.4826, a constant linked to the assumption of normality of the data, disregarding the abnormality induced by outliers (Rousseeuw & Croux, 1993).
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
#' outliers_mad(x, threshold = 3.0, replace_outlier_value = NA,
#' show_mad_values = FALSE, show_outlier_indices = FALSE,
#' b_constant = 1.4826, digits = 2, debug = FALSE)
#'
#' @examples
#' x <- c(1, 3, 3, 6, 8, 10, 10, 1000, -1000) # 1000 is an outlier
#' outliers_mad(x)
#' outliers_mad(x, threshold = 3.0)
#' outliers_mad(x, threshold = 2.5, replace_outlier_value = -999)
#' outliers_mad(x, threshold = 1.5, show_outlier_indices = TRUE)
#' outliers_mad(x, threshold = 1.5, show_mad_values = TRUE)
#' outliers_mad(x, threshold = 1.5, show_mad_values = TRUE, replace_outlier_value = -88)
outliers_mad <- function(x, threshold = 3.0, replace_outlier_value = NA, show_mad_values = FALSE, show_outlier_indices = FALSE, b_constant = 1.4826, digits = 2, debug = FALSE) {

  mad_deviate <- (x - stats::median(x, na.rm = T)) / stats::mad(x, constant = b_constant, na.rm = T)
  abs_mad_deviate <- abs(mad_deviate)

  x[abs_mad_deviate > threshold] <- replace_outlier_value
  outlier_indices <- which(abs_mad_deviate > threshold)

  if (debug) message(paste0(length(outlier_indices), " outliers detected."))

  if (show_mad_values & show_outlier_indices) {
    warning("When show_mad_values and show_outliers_indices are TRUE, returns MAD values.")
    return(round(mad_deviate, digits))
  }

  if (show_mad_values) {
    return(round(mad_deviate, digits))
  }

  if (show_outlier_indices) {
    return(outlier_indices)
  }

  return(round(x, digits))
}
