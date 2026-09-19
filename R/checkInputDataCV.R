#' @title Check Input Data for cv.srlars Function (Internal)
#'
#' @description
#' Internal helper function to validate the arguments specific to \code{cv.srlars}
#' (\code{share_grid} and \code{outer_folds}), on top of the shared \code{srlars} argument
#' checks performed by \code{checkInputData}.
#'
#' @param share_grid Integer vector of candidate \code{max_share} values.
#' @param n_models Number of models in the ensemble.
#' @param outer_folds Number of outer cross-validation folds.
#'
#' @return NULL. Stops execution with an error message if invalid inputs are detected.
#'
#' @keywords internal
#'
checkInputDataCV <- function(share_grid, n_models, outer_folds) {

    # 1. Checking share_grid
    if (!inherits(share_grid, "numeric") && !inherits(share_grid, "integer")) {
        stop("share_grid should be numeric")
    }
    if (any(share_grid != floor(share_grid))) {
        stop("share_grid should contain integers")
    }
    if (any(share_grid < 1) || any(share_grid > n_models)) {
        stop("share_grid values should be between 1 and n_models")
    }
    if (any(duplicated(share_grid))) {
        stop("share_grid should not contain duplicate values")
    }

    # 2. Checking outer_folds
    if (!inherits(outer_folds, "numeric") || outer_folds <= 1 || outer_folds != floor(outer_folds)) {
        stop("outer_folds should be a positive integer greater than 1")
    }
}
