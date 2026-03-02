#' FSCRE Competitive Selection Loop (Internal)
#'
#' @description
#' Implements Stage 2 of the FSCRE algorithm: the iterative competitive selection.
#'
#' @param Rx Global robust predictor correlation matrix.
#' @param ry Global robust predictor-response correlation vector.
#' @param x.imp Imputed design matrix (for CV).
#' @param y.imp Imputed response vector (for CV).
#' @param n_models Number of models (K).
#' @param max_predictors Maximum total predictors to select.
#' @param tolerance Stopping tolerance.
#'
#' @return A list containing the final active sets.
#'
#' @keywords internal
#'
performSelectionLoop <- function(Rx, ry,
                                 x.imp, y.imp,
                                 n_models,
                                 max_predictors,
                                 tolerance) {

    p <- ncol(x.imp)

    # --- Initialization ---

    active.sets <- vector("list", n_models)
    sign.vectors <- vector("list", n_models)
    current.correlations <- vector("list", n_models)

    for(k in 1:n_models) {
        active.sets[[k]] <- integer(0)
        sign.vectors[[k]] <- integer(0)
        current.correlations[[k]] <- ry
    }

    available.vars <- 1:p

    current.cv.errors <- numeric(n_models)
    empty_error <- computeCVError(x.imp, y.imp, integer(0))
    for(k in 1:n_models) current.cv.errors[k] <- empty_error

    n.selected <- 0
    continue.selection <- TRUE

    # --- Main Loop ---

    while (continue.selection && n.selected < max_predictors && length(available.vars) > 0) {

        # 1. Propose Candidates
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

        # 2. Evaluate Benefits
        benefits <- numeric(n_models)
        benefits[] <- -Inf

        for (k in 1:n_models) {
            cand <- candidates[[k]]
            if (!is.null(cand$next_var)) {
                new.error <- computeCVError(x.imp, y.imp, c(active.sets[[k]], cand$next_var))
                benefits[k] <- current.cv.errors[k] - new.error
            }
        }

        # Check if any valid proposal exists
        if (all(benefits == -Inf)) {
            break
        }

        # 3. Find Winner
        max.ben <- max(benefits)

        # Random Tie-Breaking
        tied.models <- which(benefits == max.ben)
        if (length(tied.models) > 1) {
            winner.k <- sample(tied.models, 1)
        } else {
            winner.k <- tied.models
        }

        winner.cand <- candidates[[winner.k]]
        winner.var <- winner.cand$next_var
        winner.base.error <- current.cv.errors[winner.k]

        # 4. Check Stopping Criterion
        ratio <- max.ben / winner.base.error

        if (ratio > tolerance) {

            # 5. Update State
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

    return(list(active.sets = active.sets))
}
