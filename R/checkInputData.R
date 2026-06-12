#' @title Check Input Data for srlars Function
#'
#' @description
#' Internal helper function to validate arguments passed to \code{srlars}.
#' Checks types, dimensions, and logical constraints.
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param n_models Number of models in the ensemble.
#' @param tolerance Relative improvement tolerance for stopping.
#' @param max_predictors Maximum total number of variables to select.
#' @param x_preprocess X cleaning method.
#' @param y_preprocess y cleaning method.
#' @param cor_estimator Correlation method.
#' @param cv_preprocess Foldwise or global CV.
#' @param cv_loss Arbiter loss function.
#' @param cv_fit Arbiter fit function.
#' @param cv_folds Number of CV folds.
#' @param compute_coef Logical.
#'
#' @return NULL. Stops execution with an error message if invalid inputs are detected.
#'
#' @keywords internal
#'
checkInputData <- function(x, y,
                           n_models,
                           tolerance,
                           max_predictors,
                           x_preprocess,
                           y_preprocess,
                           cor_estimator,
                           cv_preprocess,
                           cv_loss,
                           cv_fit,
                           cv_folds,
                           compute_coef) {

    # 1. Checking x and y
    if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
        stop("x should belong to one of the following classes: matrix, data.frame")
    }

    if (all(!inherits(y, "matrix"), !inherits(y, "numeric"))) {
        stop("y should belong to one of the following classes: matrix, numeric")
    }

    # If NO preprocessing is selected, check for NAs/NaNs/Inf.
    # If preprocessing IS selected, allow them because DDC/wrap can handle them.
    if (x_preprocess == "none") {
        if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
            stop("x should not have missing, infinite or nan values when x_preprocess='none'")
        }
    }
    
    if (y_preprocess == "none") {
        if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
            stop("y should not have missing, infinite or nan values when y_preprocess='none'")
        }
    }

    # Ensure y is a vector and dimensions match
    if(inherits(y, "matrix")) {
        if (ncol(y) > 1){
            stop("y should be a vector")
        }
    }

    if (length(y) != nrow(x)) {
        stop("y and x should have the same number of rows")
    }

    # 2. Checking n_models (K)
    if (!is.null(n_models)) {
        if (!inherits(n_models, "numeric")) {
            stop("n_models should be numeric")
        } else if (any(!n_models == floor(n_models), n_models <= 0)) {
            stop("n_models should be a positive integer greater than 0")
        }
    }

    # 3. Checking max_predictors
    if (!is.null(max_predictors)) {
        if (!inherits(max_predictors, "numeric")) {
            stop("max_predictors should be numeric")
        } else if (any(!max_predictors == floor(max_predictors), max_predictors <= 0)) {
            stop("max_predictors should be a positive integer greater than 0")
        }
        if (max_predictors > ncol(x)) {
            warning("max_predictors is larger than the number of variables in x. It will be capped at p.")
        }
    }
    
    # 4. Checking cv_folds
    if (!inherits(cv_folds, "numeric") || cv_folds <= 1 || cv_folds != floor(cv_folds)) {
        stop("cv_folds should be a positive integer greater than 1")
    }

    # 5. Check logical flags
    if (!(compute_coef %in% c(TRUE, FALSE)))
        stop("compute_coef should be TRUE or FALSE.")
}