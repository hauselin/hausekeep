#' @title Specification curve analysis (SCA) for linear regression
#' @name sca_lm
#'
#' @description sca_lm is used to run specification curve analysis using linear regression. Fits every possible model for your specified dependent variable, predictor variables, and covariates
#'
#' @param data dataframe
#' @param dv dependent variable (outcome variable) (character)
#' @param ivs independent variable(s) or predictor(s) (character)
#' @param covariates covariates (character)
#'
#' @return A tibble/data.table with results
#'
#' @details
#' Details to follow

#' @references
#' References to follow
#' \cr \cr
#'
#' @author Hause Lin
#'
#' @usage
#' sca_lm(data, dv, ivs, covariates = NULL)
#'
#' @importFrom dplyr left_join distinct tibble bind_rows
#' @importFrom combinat permn
#' @importFrom data.table setDT
#' @importFrom stats lm as.formula
#' @export
#' @examples
#' # model with 1 covariate
#' m1 <- sca_lm(data = mtcars, dv = "mpg", ivs = c("cyl", "carb"), covariates = c("vs"))
#' m1[, c("modelformula", "term", "estimate", "p.value")] # model, term, beta, p value
#'
#' # model without covariates
#' m2 <- sca_lm(data = mtcars, dv = "mpg", ivs = c("cyl", "gear"))
#' m2[, c("modelformula", "term", "estimate", "es.r")] # es.r is effect size
sca_lm <- function(data, dv, ivs, covariates = NULL) {

  # ivs
  n_ivs <- length(ivs)
  mat <- matrix(rep(1, n_ivs), ncol = n_ivs)
  for (i in 1:n_ivs) {
    tempv <- c(rep(1, i), rep(0, n_ivs - i))
    mattemp <- matrix(unlist(permn(tempv)), ncol = n_ivs)
    mat <- rbind(mat, mattemp)
  }
  colnames(mat) <- paste0(ivs)
  mat_iv_df <- dplyr::distinct(data.table(mat))
  mat_iv_df

  # covariates
  if (!is.null(covariates)) {
    n_covariates <- length(covariates)
    mat <- matrix(c(rep(1, n_covariates), rep(0, n_covariates)), ncol = n_covariates, byrow = T)
    for (i in 1:n_covariates) {
      tempv <- c(rep(1, i), rep(0, n_covariates - i))
      mattemp <- matrix(unlist(combinat::permn(tempv)), ncol = n_covariates)
      mat <- rbind(mat, mattemp)
    }
    colnames(mat) <- paste0(covariates)
    mat_cov_df <- dplyr::distinct(data.table(mat))
    mat_cov_df
  }

  # create all combinations of ivs and covaraites
  if (!is.null(covariates)) {
    dt1 <- tibble()
    for (i in 1:nrow(mat_iv_df)) {
      dt1 <- dplyr::bind_rows(dt1, cbind(mat_iv_df[i, ], mat_cov_df))
    }
  } else {
    dt1 <- data.table::copy(mat_iv_df)
  }
  data.table::setDT(dt1)

  # build model formula
  dt1[, modelformula := ""]
  for (i in 1:nrow(dt1)) {
    dt1[i, modelformula := paste0(dv, " ~ ", paste0(names(dt1)[which(dt1[i, ] == 1)], collapse = " + "))]
  }

  # fit models
  results <- dt1[, summaryh(stats::lm(stats::as.formula(modelformula), data = data), showTable = T)$results2,
                 by = modelformula]

  results <- dplyr::left_join(dt1, results, by = "modelformula")
  setDT(results)
  results
}
