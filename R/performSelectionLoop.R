#' @title FSCRE Competitive Selection Loop (Internal)
#'
#' @description
#' Implements Stage 2 of the FSCRE algorithm: the iterative competitive selection
#' using a proposer–arbiter mechanism. Candidate moves are proposed by a robust
#' LARS step computed from the (robust) correlation inputs, and accepted by an
#' arbiter based on cross-validated predictive improvement.
#'
#' The function supports two CV preprocessing modes:
#' \itemize{
#'   \item \code{cv_preprocess = "global"}: CV is computed on the globally preprocessed
#'   training data \code{(x.imp, y.imp)}.
#'   \item \code{cv_preprocess = "foldwise"}: preprocessing is fitted on each CV fold's
#'   training subset and applied to that fold's validation subset (via
#'   \code{cellWise::DDCpredict} for \code{x_preprocess="ddc"} and by reusing the training
#'   fold's location/scale for \code{y_preprocess="wrap"}).
#' }
#'
#' The selection step terminates when no candidate yields a strictly positive CV
#' improvement, or when the best relative improvement falls below \code{tolerance}.
#'
#' @param Rx Global predictor correlation matrix used by the LARS proposer (p x p).
#' @param ry Global predictor-response correlation vector used by the LARS proposer (length p).
#' @param x Raw design matrix (n x p). Used for foldwise preprocessing if requested.
#' @param y Raw response vector (length n). Used for foldwise preprocessing if requested.
#' @param x.imp Globally preprocessed design matrix (n x p), used when \code{cv_preprocess="global"}.
#' @param y.imp Globally preprocessed response vector (length n), used when \code{cv_preprocess="global"}.
#' @param n_models Number of models in the ensemble (K).
#' @param max_predictors Maximum total number of predictors to select across all models.
#' @param tolerance Relative improvement tolerance for stopping (\eqn{\tau}).
#' @param x_preprocess Character. Preprocessing method for predictors in the CV loop
#'   (e.g., \code{"ddc"} or \code{"none"}).
#' @param y_preprocess Character. Preprocessing method for response in the CV loop
#'   (e.g., \code{"wrap"}, \code{"robust_z"}, or \code{"none"}).
#' @param cv_preprocess Character. CV preprocessing mode: \code{"global"} or \code{"foldwise"}.
#' @param cv_fit Character. Inner CV fitting method used by \code{computeCVError}
#'   (e.g., \code{"ls"} or \code{"huber"}).
#' @param cv_loss Character. CV scoring loss used by \code{computeCVError}
#'   (e.g., \code{"mse"}, \code{"trimmed"}, or \code{"huber"}).
#' @param cv_folds Integer. Number of CV folds.
#'
#' @return A list with component:
#' \describe{
#'   \item{\code{active.sets}}{A list of length \code{n_models}, where each element is an integer
#'   vector of selected variable indices for that sub-model.}
#' }
#'
#' @keywords internal
#'
#' @importFrom stats median mad sd
#' @importFrom cellWise DDC DDCpredict wrap
#'
#' @seealso \code{\link{getLarsProposal}}, \code{\link{computeCVError}}, \code{\link{computeRobustFoundation}}
#' 
performSelectionLoop <- function(Rx, ry,
                                 x, y, x.imp, y.imp,
                                 n_models,
                                 max_predictors,
                                 tolerance,
                                 x_preprocess,
                                 y_preprocess,
                                 cv_preprocess,
                                 cv_fit,
                                 cv_loss,
                                 cv_folds) {

    n <- nrow(x)
    p <- ncol(x)

    # Ensure colnames exist for DDC/DDCpredict consistency
    if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(p))

    # ____________________________
    # 0. Setup and Cache CV Folds
    # ____________________________

    fold_ids <- sample(rep(1:cv_folds, length.out = n))
    cv_data <- vector("list", cv_folds)

    for (f in 1:cv_folds) {
        train_idx <- which(fold_ids != f)
        val_idx   <- which(fold_ids == f)

        if (cv_preprocess == "foldwise") {

            # ___________________________________________________________
            # (A) X preprocessing: fit on fold-train, apply to fold-val
            # ___________________________________________________________
          
            x_train_raw <- x[train_idx, , drop = FALSE]
            x_val_raw   <- x[val_idx,   , drop = FALSE]

            if (x_preprocess == "ddc") {
                ddc_train <- cellWise::DDC(
                    x_train_raw,
                    DDCpars = list(fastDDC = TRUE, silent = TRUE)
                )
                x_train_f <- ddc_train$Ximp

                # Apply the fitted DDC to validation fold
                ddc_val <- cellWise::DDCpredict(x_val_raw, ddc_train)
                x_val_f <- ddc_val$Ximp
            } else {
                x_train_f <- x_train_raw
                x_val_f   <- x_val_raw
            }

            # ___________________________________________________________
            # (B) y preprocessing: fit on fold-train, apply to fold-val
            # ___________________________________________________________
          
            y_train_raw <- as.numeric(y[train_idx])
            y_val_raw   <- as.numeric(y[val_idx])

            if (y_preprocess == "wrap") {
                w_train <- cellWise::wrap(as.matrix(y_train_raw))
                y_train_f <- as.numeric(w_train$Xw)

                # Apply same loc/scale to val fold to keep scales consistent
                w_val <- cellWise::wrap(as.matrix(y_val_raw),
                                        locX = w_train$loc,
                                        scaleX = w_train$scale)
                y_val_f <- as.numeric(w_val$Xw)

            } else if (y_preprocess == "robust_z") {
                med_y <- stats::median(y_train_raw)
                mad_y <- stats::mad(y_train_raw)
                if (mad_y < 1e-12) mad_y <- stats::sd(y_train_raw)
                if (mad_y < 1e-12) mad_y <- 1

                clip_lo <- med_y - 3 * mad_y
                clip_hi <- med_y + 3 * mad_y

                y_train_f <- pmin(pmax(y_train_raw, clip_lo), clip_hi)
                y_val_f   <- pmin(pmax(y_val_raw,   clip_lo), clip_hi)

            } else {
                y_train_f <- y_train_raw
                y_val_f   <- y_val_raw
            }

        } else {
            # Global: split already-preprocessed data (fast, but not strictly cross-fitted)
            x_train_f <- x.imp[train_idx, , drop = FALSE]
            y_train_f <- y.imp[train_idx]
            x_val_f   <- x.imp[val_idx,   , drop = FALSE]
            y_val_f   <- y.imp[val_idx]
        }

        cv_data[[f]] <- list(
            x_train = x_train_f,
            y_train = as.numeric(y_train_f),
            x_val   = x_val_f,
            y_val   = as.numeric(y_val_f)
        )
    }

    # ______________________
    # 1. Initialization
    # ______________________

    active.sets <- vector("list", n_models)
    sign.vectors <- vector("list", n_models)
    current.correlations <- vector("list", n_models)

    for (k in 1:n_models) {
        active.sets[[k]] <- integer(0)
        sign.vectors[[k]] <- integer(0)
        current.correlations[[k]] <- ry
    }

    available.vars <- 1:p

    current.cv.errors <- numeric(n_models)
    empty_error <- computeCVError(cv_data, integer(0), cv_fit, cv_loss)
    for (k in 1:n_models) current.cv.errors[k] <- empty_error

    n.selected <- 0
    continue.selection <- TRUE

    # ______________________________
    # 2. Main Proposer-Arbiter Loop
    # ______________________________

    while (continue.selection && n.selected < max_predictors && length(available.vars) > 0) {

        # A. Propose Candidates
        candidates <- vector("list", n_models)
        for (k in 1:n_models) {
            candidates[[k]] <- getLarsProposal(
                Rx,
                active.sets[[k]],
                sign.vectors[[k]],
                current.correlations[[k]],
                available.vars
            )
        }

        # B. Evaluate Benefits
        benefits <- rep(-Inf, n_models)
        for (k in 1:n_models) {
            cand <- candidates[[k]]
            if (!is.null(cand$next_var)) {
                new.error <- computeCVError(cv_data,
                                            c(active.sets[[k]], cand$next_var),
                                            cv_fit, cv_loss)
                benefits[k] <- current.cv.errors[k] - new.error
            }
        }

        if (all(benefits == -Inf)) break

        # C. Winner
        max.ben <- max(benefits)
        tied.models <- which(benefits == max.ben)
        winner.k <- if (length(tied.models) > 1) sample(tied.models, 1) else tied.models

        winner.cand <- candidates[[winner.k]]
        winner.var  <- winner.cand$next_var
        winner.base.error <- current.cv.errors[winner.k]

        # --- NEW: enforce positive benefit (Algorithm 1: B* > 0) ---
        if (!is.finite(max.ben) || max.ben <= 0) {
            break
        }

        # --- NEW: guard against invalid baseline error ---
        if (!is.finite(winner.base.error) || winner.base.error <= 0) {
            break
        }

        ratio <- max.ben / winner.base.error

        if (ratio > tolerance) {
            active.sets[[winner.k]] <- c(active.sets[[winner.k]], winner.var)
            sign.vectors[[winner.k]] <- c(sign.vectors[[winner.k]], winner.cand$next_sign)

            current.cv.errors[winner.k] <- winner.base.error - max.ben
            current.correlations[[winner.k]] <- current.correlations[[winner.k]] -
                (winner.cand$gamma * winner.cand$a_vec)

            available.vars <- setdiff(available.vars, winner.var)
            n.selected <- n.selected + 1
        } else {
            continue.selection <- FALSE
        }
    }

    list(active.sets = active.sets)
}