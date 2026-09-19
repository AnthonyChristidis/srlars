#' @title Ensemble Prediction from a Raw Fit List (Internal)
#'
#' @description
#' Averages the per-model linear predictions of an \code{srlars}-style fit list
#' (\code{coefficients}, \code{intercepts}) over a design matrix \code{X}. Used by
#' \code{\link{cv.srlars}} to score candidate \code{max_share} values without going through
#' \code{predict.srlars} (which expects a full \code{srlars} object, not an intermediate
#' candidate fit).
#'
#' @param X Design matrix.
#' @param fit A list with components \code{coefficients} (list of numeric vectors) and
#'   \code{intercepts} (numeric vector), as returned by \code{computeFinalFit}.
#'
#' @return Numeric vector of ensemble-averaged predictions.
#'
#' @keywords internal
#'
srlars_ensemble_predict <- function(X, fit) {

    n_models <- length(fit$coefficients)
    preds <- numeric(nrow(X))

    for (k in 1:n_models) {
        preds <- preds + (fit$intercepts[k] + as.numeric(X %*% fit$coefficients[[k]]))
    }

    preds / n_models
}

#' @importFrom stats coef sd median mad cor predict rnorm runif rbinom
#' @importFrom cellWise DDC wrap DDCpredict
#' @importFrom robustbase lmrob
#'
#' @title Cross-Validated Selection of \code{max_share} for FSCRE (cv.srlars)
#'
#' @description
#' \code{cv.srlars} chooses the \code{max_share} argument of \code{\link{srlars}} by
#' cross-validation and returns the FSCRE ensemble refit on the full data at the
#' cross-validated optimum. The internal cross-validation used by the competitive arbiter
#' inside \code{\link{srlars}} (controlled by \code{cv_folds}) is unrelated and untouched by
#' this function: it decides which variable each sub-model proposes next. \code{cv.srlars}
#' adds a separate, outer cross-validation loop whose only job is to pick the best
#' \code{max_share} value from \code{share_grid}.
#'
#' To keep this efficient, the expensive cellwise-robust preprocessing stage
#' (DDC imputation, wrapping, and robust correlation estimation) does not depend on
#' \code{max_share}, so it is computed only once per outer fold and reused across every
#' candidate value in \code{share_grid} for that fold, rather than being recomputed for every
#' (fold, candidate) pair.
#'
#' @param x Design matrix (n x p).
#' @param y Response vector (n x 1).
#' @param n_models Number of models in the ensemble (K). Default is 5.
#' @param share_grid Integer vector of candidate \code{max_share} values to evaluate. Default is
#' \code{1:n_models}.
#' @param outer_folds Integer. Number of outer cross-validation folds used to select
#' \code{max_share}. Default is 5.
#' @param tolerance Relative improvement tolerance for stopping (tau). Default is 1e-8.
#' @param n_min Integer or NULL. Minimum number of variables each sub-model is guaranteed
#' (subject to availability), passed through unchanged to every candidate fit and the final
#' refit -- \code{cv.srlars} does not tune \code{n_min}, only \code{max_share}. Default is NULL
#' (no floor). See \code{\link{srlars}}.
#' @param max_predictors Maximum total number of variables to select across all models. Default is n * n_models.
#' @param x_preprocess Character. "ddc" (default) for cellwise cleaning, or "none".
#' @param y_preprocess Character. "wrap" (default) for univariate robustification, "robust_z", or "none".
#' @param cor_estimator Character. "wrap" (default) for robust PSD correlation, or "pearson".
#' @param cv_preprocess Character. "global" (default) or "foldwise" (to prevent data leakage).
#' @param cv_fit Character. "huber" (default) or "ls" for the inner arbiter fitting method.
#' @param cv_loss Character. "huber" (default), "trimmed", or "mse" for arbiter scoring and for
#' scoring \code{max_share} candidates in the outer loop.
#' @param cv_folds Integer. Number of internal (arbiter) cross-validation folds. Default is 5.
#' @param compute_coef Logical. If TRUE, fits the final robust MM-models. Default is TRUE.
#'
#' @return An object of class \code{c("cv.srlars", "srlars")}: the \code{srlars} fit at the
#' cross-validated optimal \code{max_share}, with three extra components:
#' \describe{
#'   \item{\code{max_share}}{The cross-validated optimal value.}
#'   \item{\code{share_grid}}{The candidate values that were evaluated.}
#'   \item{\code{cv_errors}}{The mean out-of-fold error for each value in \code{share_grid}.}
#' }
#' Because no \code{coef.cv.srlars} or \code{predict.cv.srlars} methods are defined, calling
#' \code{coef()} or \code{predict()} on the result dispatches to \code{\link{coef.srlars}} /
#' \code{\link{predict.srlars}} and operates on this optimal fit directly.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @seealso \code{\link{srlars}}, \code{\link{coef.srlars}}, \code{\link{predict.srlars}}
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
#' p <- 30
#' rho.within <- 0.8
#' rho.between <- 0.2
#' p.active <- 10
#' group.size <- 5
#' snr <- 3
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
#' colnames(x) <- paste0("V", 1:p)
#' y <- x %*% true.beta + rnorm(n, 0, sigma)
#'
#' # Cross-validated choice of max_share
#' cv_fit <- cv.srlars(x, y,
#'                     n_models = 3,
#'                     outer_folds = 3,
#'                     tolerance = 1e-4)
#'
#' # Cross-validated optimal max_share
#' print(cv_fit$max_share)
#'
#' # coef() and predict() dispatch to the srlars methods automatically
#' cv_coefs <- coef(cv_fit)
#'
cv.srlars <- function(x, y,
                      n_models = 5,
                      share_grid = NULL,
                      outer_folds = 5,
                      tolerance = 1e-8,
                      n_min = NULL,
                      max_predictors = NULL,
                      x_preprocess = c("ddc", "none"),
                      y_preprocess = c("wrap", "robust_z", "none"),
                      cor_estimator = c("wrap", "pearson"),
                      cv_preprocess = c("global", "foldwise"),
                      cv_fit = c("huber", "ls"),
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

    if (is.null(share_grid)) {
        share_grid <- 1:n_models
    }

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
                   compute_coef,
                   max_share = NULL,
                   n_min = n_min)
    checkInputDataCV(share_grid, n_models, outer_folds)

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

    if (is.null(colnames(x))) {
        colnames(x) <- paste0("V", 1:ncol(x))
    }

    # ________________________________________________________
    # 3. Outer CV Loop: Score Each max_share Candidate per Fold
    # ________________________________________________________

    outer.fold.ids <- sample(rep(1:outer_folds, length.out = n))
    grid.errors <- matrix(NA_real_, nrow = outer_folds, ncol = length(share_grid))

    for (f in 1:outer_folds) {

        train.idx <- which(outer.fold.ids != f)
        val.idx   <- which(outer.fold.ids == f)

        x.train <- x[train.idx, , drop = FALSE]
        y.train <- y[train.idx]
        x.val   <- x[val.idx,   , drop = FALSE]
        y.val   <- y[val.idx] # Raw, unprocessed response used for scoring

        # Robust foundation computed ONCE per outer fold, reused across the whole share_grid
        foundation <- computeRobustFoundation(x.train, y.train, x_preprocess, y_preprocess, cor_estimator)

        # Clean the held-out predictors into the same feature space using this fold's own DDC fit
        x.val.clean <- x.val
        if (x_preprocess == "ddc" && !is.null(foundation$ddc.object)) {
            x.val.clean <- tryCatch({
                aug <- cbind(x.val, NA)
                colnames(aug) <- c(colnames(x.val), ".y_placeholder")
                ddc.pred <- cellWise::DDCpredict(Xnew = aug, InitialDDC = foundation$ddc.object)
                ddc.pred$Ximp[, seq_len(ncol(x.val)), drop = FALSE]
            }, error = function(e) x.val)
        }

        for (g in seq_along(share_grid)) {

            share.g <- share_grid[g]

            selection.g <- performSelectionLoop(foundation$Rx, foundation$ry,
                                                x.train, y.train, foundation$x.imp, foundation$y.imp,
                                                n_models, max_predictors, tolerance,
                                                x_preprocess, y_preprocess,
                                                cv_preprocess, cv_fit, cv_loss, cv_folds,
                                                max_share = share.g, n_min = n_min)

            fit.g <- computeFinalFit(foundation$x.imp, foundation$y.imp, selection.g$active.sets,
                                     compute_coef = TRUE)

            train.pred <- srlars_ensemble_predict(foundation$x.imp, fit.g)
            val.pred   <- srlars_ensemble_predict(x.val.clean, fit.g)

            r.train <- foundation$y.imp - train.pred
            r.val   <- y.val - val.pred

            grid.errors[f, g] <- computeRobustLoss(r.val, if (cv_loss == "huber") r.train else NULL, cv_loss)
        }
    }

    mean.errors <- colMeans(grid.errors)
    best.share <- share_grid[which.min(mean.errors)]

    # ______________________________________________
    # 4. Final Refit on Full Data at the CV-Optimum
    # ______________________________________________

    final.fit <- srlars(x, y,
                        n_models = n_models,
                        max_share = best.share,
                        tolerance = tolerance,
                        n_min = n_min,
                        max_predictors = max_predictors,
                        x_preprocess = x_preprocess,
                        y_preprocess = y_preprocess,
                        cor_estimator = cor_estimator,
                        cv_preprocess = cv_preprocess,
                        cv_fit = cv_fit,
                        cv_loss = cv_loss,
                        cv_folds = cv_folds,
                        compute_coef = compute_coef)

    final.fit$max_share <- best.share
    final.fit$share_grid <- share_grid
    final.fit$cv_errors <- mean.errors

    class(final.fit) <- c("cv.srlars", class(final.fit))
    return(final.fit)
}
