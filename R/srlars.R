#' @importFrom stats coef sd median mad cor predict rnorm runif rbinom
#' @importFrom cellWise DDC wrap
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
#' @param x_preprocess Character. "ddc" (default) for cellwise cleaning, or "none".
#' @param y_preprocess Character. "wrap" (default) for univariate robustification, "robust_z", or "none".
#' @param cor_estimator Character. "wrap" (default) for robust PSD correlation, or "pearson".
#' @param cv_preprocess Character. "global" (default) or "foldwise" (to prevent data leakage).
#' @param cv_fit Character. "ls" (default) or "huber" for the inner arbiter fitting method.
#' @param cv_loss Character. "huber" (default), "trimmed", or "mse" for arbiter scoring.
#' @param cv_folds Integer. Number of cross-validation folds. Default is 5.
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
srlars <- function(x, y,
                   n_models = 5,
                   tolerance = 1e-8,
                   max_predictors = NULL,
                   x_preprocess = c("ddc", "none"),
                   y_preprocess = c("wrap", "robust_z", "none"),
                   cor_estimator = c("wrap", "pearson"),
                   cv_preprocess = c("global", "foldwise"),
                   cv_fit = c("ls", "huber"),
                   cv_loss = c("huber", "trimmed", "mse"),
                   cv_folds = 5,
                   compute_coef = TRUE) {

    # Match arguments to ensure valid inputs
    x_preprocess <- match.arg(x_preprocess)
    y_preprocess <- match.arg(y_preprocess)
    cor_estimator <- match.arg(cor_estimator)
    cv_preprocess <- match.arg(cv_preprocess)
    cv_fit <- match.arg(cv_fit)
    cv_loss <- match.arg(cv_loss)

    # ________________
    # 1. Input Checks
    # ________________
  
    checkInputData(x, y,
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
                   compute_coef)
    
    # _________
    # 2. Setup
    # _________
  
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

    # ______________________________________________________________________________
    # 3. Stage 1: Global Robust Foundation (For LARS initialization and Final Fit)
    # ______________________________________________________________________________
  
    # X and Y are now processed separately to avoid target leakage.
    foundation <- computeRobustFoundation(x, y, x_preprocess, y_preprocess, cor_estimator)

    x.imp <- foundation$x.imp
    y.imp <- foundation$y.imp
    Rx <- foundation$Rx
    ry <- foundation$ry
    ddc.object <- foundation$ddc.object

    # ________________________________________
    # 4. Stage 2: Competitive Selection Loop
    # ________________________________________
  
    # Passes all new CV logic down to the arbiter.
    selection.results <- performSelectionLoop(Rx, ry, x, y, x.imp, y.imp,
                                              n_models, max_predictors, tolerance,
                                              x_preprocess, y_preprocess,
                                              cv_preprocess, cv_fit, cv_loss, cv_folds)

    # ______________________
    # 5. Stage 3: Final Fit
    # ______________________
  
    # Fit MM-estimators on the globally cleaned data
    final.model <- computeFinalFit(x.imp, y.imp, selection.results$active.sets, compute_coef)

    # ________________________
    # 6. Output Construction
    # ________________________
  
    object <- list(
        active.sets = selection.results$active.sets,
        coefficients = final.model$coefficients,
        intercepts = final.model$intercepts,
        n_models = n_models,
        x_preprocess = x_preprocess,
        y_preprocess = y_preprocess,
        ddc.object = ddc.object,
        imputed.data = list(x = x.imp, y = y.imp) # Always return this for user transparency
    )

    class(object) <- "srlars"
    return(object)
}