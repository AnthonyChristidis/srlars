#' @importFrom stats coef sd median mad cor predict rnorm runif rbinom
#' @importFrom cellWise DDC
#' @importFrom robustbase lmrob
#' @importFrom mvnfast rmvn
#'
#' @title Fast and Scalable Cellwise-Robust Ensemble (FSCRE)
#'
#' @description
#' \code{srlars} performs the FSCRE algorithm for robust variable selection and regression.
#'
#' @param x Design matrix (n x p).
#' @param y Response vector (n x 1).
#' @param n_models Number of models in the ensemble (K). Default is 5.
#' @param tolerance Relative improvement tolerance for stopping (tau). Default is 1e-8.
#' @param max_predictors Maximum total number of variables to select across all models. Default is n * n_models.
#' @param robust Logical. If TRUE (default), performs DDC imputation and robust initialization.
#' @param compute_coef Logical. If TRUE, fits the final robust MM-models. Default is TRUE.
#'
#' @return An object of class \code{srlars} containing the selected variables and coefficients.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @seealso \code{\link{coef.srlars}}, \code{\link{predict.srlars}}
#'
#' @export
#'
#' @examples
#' # Required libraries
#' library(mvnfast)
#' library(cellWise)
#' library(robustbase)
#'
#' # Simulation parameters
#' n <- 50
#' p <- 100
#' rho.within <- 0.8
#' rho.between <- 0.2
#' p.active <- 20
#' group.size <- 5
#' snr <- 3
#' contamination.prop <- 0.1
#'
#' # Setting the seed
#' set.seed(0)
#'
#' # Block correlation structure
#' sigma.mat <- matrix(0, p, p)
#' sigma.mat[1:p.active, 1:p.active] <- rho.between
#' for(group in 0:(p.active/group.size - 1))
#'   sigma.mat[(group*group.size+1):(group*group.size+group.size),
#'   (group*group.size+1):(group*group.size+group.size)] <- rho.within
#' diag(sigma.mat) <- 1
#'
#' # Simulation of beta vector
#' true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
#'
#' # Setting the SD of the variance
#' sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))
#'
#' # Simulation of uncontaminated data
#' x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
#' y <- x %*% true.beta + rnorm(n, 0, sigma)
#'
#' # Cellwise contamination
#' contamination_indices <- sample(1:(n * p), round(n * p * contamination.prop))
#' x_train <- x
#' x_train[contamination_indices] <- runif(length(contamination_indices), -10, 10)
#'
#' # FSCRE Ensemble model
#' ensemble_fit <- srlars(x_train, y,
#'                        n_models = 5,
#'                        tolerance = 1e-8,
#'                        robust = TRUE,
#'                        compute_coef = TRUE)
#'
#' # Check selected variables
#' print(ensemble_fit$active.sets)
#'
srlars <- function(x, y,
                   n_models = 5,
                   tolerance = 1e-8,
                   max_predictors = NULL,
                   robust = TRUE,
                   compute_coef = TRUE) {

    # 1. Input Checks
    checkInputData(x, y,
                   n_models,
                   tolerance,
                   max_predictors,
                   robust,
                   compute_coef)

    # 2. Setup
    n <- nrow(x)
    p <- ncol(x)
    x <- as.matrix(x)
    y <- as.numeric(y)

    if (is.null(max_predictors)) {
        max_predictors <- min(p, n * n_models)
    }

    # Ensure Column Names (Critical for DDC matching in predict)
    if (is.null(colnames(x))) {
        colnames(x) <- paste0("V", 1:ncol(x))
    }

    # 3. Stage 1: Robust Foundation
    # We use the internal helper to perform DDC and compute correlations
    foundation <- computeRobustFoundation(x, y, robust)

    x.imp <- foundation$x.imp
    y.imp <- foundation$y.imp
    Rx <- foundation$Rx
    ry <- foundation$ry

    # Capture the DDC object for prediction
    ddc.object <- foundation$ddc.object

    # 4. Stage 2: Competitive Selection Loop
    # The iterative "proposer-arbiter" process
    selection.results <- performSelectionLoop(Rx, ry, x.imp, y.imp,
                                              n_models, max_predictors, tolerance)

    # 5. Stage 3: Final Fit
    # Fit MM-estimators on the selected sets (if requested)
    final.model <- computeFinalFit(x.imp, y.imp, selection.results$active.sets, compute_coef)

    # 6. Output Construction
    object <- list(
        active.sets = selection.results$active.sets,
        coefficients = final.model$coefficients,
        intercepts = final.model$intercepts,
        n_models = n_models,
        robust = robust,
        ddc.object = ddc.object,
        imputed.data = if(!robust) list(x = x.imp, y = y.imp) else NULL
    )

    class(object) <- "srlars"
    return(object)
}
