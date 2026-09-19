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
#' @param max_share Integer. Maximum number of sub-models (1 to n_models) in which a given
#'   variable may appear. Default is 1 (fully disjoint sub-models). For \code{1 < max_share <
#'   n_models}, each sub-model's first selected variable is forced distinct across sub-models;
#'   sharing is only permitted afterward. This restriction is lifted when \code{max_share =
#'   n_models}.
#' @param n_min Integer or NULL. Minimum number of variables each sub-model is guaranteed
#'   (subject to availability under the \code{max_share}/diversity pool restrictions), even if no
#'   candidate clears the usual positive-benefit or \code{tolerance} requirement. Default is NULL
#'   (no floor enforced, original behavior).
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
                                 cv_folds,
                                 max_share = 1,
                                 n_min = NULL) {

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

    var.usage <- integer(p)
    seed.vars <- integer(0)

    current.cv.errors <- numeric(n_models)
    empty_error <- computeCVError(cv_data, integer(0), cv_fit, cv_loss)
    for (k in 1:n_models) current.cv.errors[k] <- empty_error

    n.selected <- 0
    continue.selection <- TRUE

    # ______________________________
    # 2. Main Proposer-Arbiter Loop
    # ______________________________

    while (continue.selection && n.selected < max_predictors && any(var.usage < max_share)) {

        # A. Propose Candidates
        candidates <- vector("list", n_models)
        for (k in 1:n_models) {
            pool.k <- which(var.usage < max_share)
            if (length(active.sets[[k]]) > 0) {
                pool.k <- setdiff(pool.k, active.sets[[k]])
            } else if (max_share < n_models) {
                # Force distinct seeds across sub-models when any sharing is allowed but not
                # unrestricted: a variable already used as another model's first pick cannot
                # be proposed as a fresh model's seed too, preventing several models from
                # redundantly duplicating the same "obviously best" cold-start variable.
                pool.k <- setdiff(pool.k, seed.vars)
            }

            candidates[[k]] <- getLarsProposal(
                Rx,
                active.sets[[k]],
                sign.vectors[[k]],
                current.correlations[[k]],
                pool.k
            )
        }

        # A2. Floor Enforcement (n_min)
        # Force-accept a candidate for a sub-model still below the floor, bypassing the
        # positive-benefit/tolerance requirement below -- but never bypassing the pool
        # restrictions already applied above (max_share usage cap and seed diversity).
        if (!is.null(n_min)) {
            below.floor <- which(vapply(active.sets, length, integer(1)) < n_min)
            below.floor <- below.floor[!vapply(candidates[below.floor],
                                               function(cand) is.null(cand$next_var),
                                               logical(1))]

            if (length(below.floor) > 0) {
                forced.benefits <- rep(-Inf, length(below.floor))
                for (i in seq_along(below.floor)) {
                    k <- below.floor[i]
                    cand <- candidates[[k]]
                    new.error <- computeCVError(cv_data,
                                                c(active.sets[[k]], cand$next_var),
                                                cv_fit, cv_loss)
                    forced.benefits[i] <- current.cv.errors[k] - new.error
                }

                if (any(is.finite(forced.benefits))) {
                    best.forced <- max(forced.benefits[is.finite(forced.benefits)])
                    best.local <- which(forced.benefits == best.forced)
                    winner.local <- if (length(best.local) > 1) sample(best.local, 1) else best.local
                    winner.k <- below.floor[winner.local]

                    winner.cand <- candidates[[winner.k]]
                    winner.var <- winner.cand$next_var
                    winner.base.error <- current.cv.errors[winner.k]
                    winner.benefit <- forced.benefits[winner.local]

                    if (length(active.sets[[winner.k]]) == 0) {
                        seed.vars <- c(seed.vars, winner.var)
                    }
                    active.sets[[winner.k]] <- c(active.sets[[winner.k]], winner.var)
                    sign.vectors[[winner.k]] <- c(sign.vectors[[winner.k]], winner.cand$next_sign)

                    current.cv.errors[winner.k] <- winner.base.error - winner.benefit
                    current.correlations[[winner.k]] <- current.correlations[[winner.k]] -
                        (winner.cand$gamma * winner.cand$a_vec)

                    var.usage[winner.var] <- var.usage[winner.var] + 1L
                    n.selected <- n.selected + 1
                    next
                }
            }
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
            if (length(active.sets[[winner.k]]) == 0) {
                seed.vars <- c(seed.vars, winner.var)
            }
            active.sets[[winner.k]] <- c(active.sets[[winner.k]], winner.var)
            sign.vectors[[winner.k]] <- c(sign.vectors[[winner.k]], winner.cand$next_sign)

            current.cv.errors[winner.k] <- winner.base.error - max.ben
            current.correlations[[winner.k]] <- current.correlations[[winner.k]] -
                (winner.cand$gamma * winner.cand$a_vec)

            var.usage[winner.var] <- var.usage[winner.var] + 1L
            n.selected <- n.selected + 1
        } else {
            continue.selection <- FALSE
        }
    }

    list(active.sets = active.sets)
}