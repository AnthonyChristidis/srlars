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
#' @param robust Logical.
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
                           robust,
                           compute_coef) {

    # 1. Checking x and y
    if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
        stop("x should belong to one of the following classes: matrix, data.frame")
    }

    if (all(!inherits(y, "matrix"), !inherits(y, "numeric"))) {
        stop("y should belong to one of the following classes: matrix, numeric")
    }

    # Note: Removed check for NAs/NaNs/Inf in x and y because
    # DDC is explicitly designed to handle them. We should allow them if robust=TRUE.
    if (!robust) {
        if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
            stop("x should not have missing, infinite or nan values when robust=FALSE")
        }
        if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
            stop("y should not have missing, infinite or nan values when robust=FALSE")
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

    # # 3. Checking tolerance (tau)
    # if (!inherits(tolerance, "numeric")) {
    #     stop("tolerance should be numeric")
    # } else if (tolerance < 0) {
    #     stop("tolerance should be non-negative")
    # }

    # 4. Checking max_predictors
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

    # 5. Check logical flags
    if (!(robust %in% c(TRUE, FALSE)))
        stop("robust should be TRUE or FALSE.")

    if (!(compute_coef %in% c(TRUE, FALSE)))
        stop("compute_coef should be TRUE or FALSE.")
}
