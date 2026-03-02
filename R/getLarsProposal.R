#' @title Get Robust LARS Proposal (Internal)
#'
#' @description
#' Calculates the next LARS step analytically using correlation matrices.
#'
#' @param Rx Global correlation matrix (p x p).
#' @param active.set Integer vector of active indices.
#' @param sign.vector Integer vector of signs for active set.
#' @param current.correlations Numeric vector of current correlations with residual.
#' @param available.vars Integer vector of candidate indices to search.
#'
#' @return A list containing: next_var, next_sign, gamma, a_vec (for update). Or NULL.
#'
#' @keywords internal
#'
getLarsProposal <- function(Rx,
                            active.set,
                            sign.vector,
                            current.correlations,
                            available.vars) {

    p <- ncol(Rx)

    # --- Case 0: Empty Active Set ---
    if (length(active.set) == 0) {
        # If empty, pick variable with max absolute correlation from available pool
        pool_corrs <- abs(current.correlations[available.vars])

        if (length(pool_corrs) == 0) return(NULL)

        max_idx <- which.max(pool_corrs)
        next_var <- available.vars[max_idx]
        next_sign <- sign(current.correlations[next_var])

        return(list(
            next_var = next_var,
            next_sign = next_sign,
            gamma = 0,       # No movement yet
            a_vec = rep(0, p) # No update vector yet
        ))
    }

    # --- Case 1: Non-Empty Active Set ---

    # 1. Get current max correlation (r_A)
    r_A <- abs(current.correlations[active.set[1]])

    # 2. Compute Geometric Quantities
    # We need (D R D)^-1.
    R_S <- Rx[active.set, active.set, drop=FALSE]

    # CRITICAL FIX: Add Ridge Regularization for Stability
    # This prevents the matrix from becoming singular due to high collinearity
    # allowing the LARS path to continue proposing variables.
    ridge_lambda <- 1e-5
    diag(R_S) <- diag(R_S) + ridge_lambda

    R_S_signed <- R_S * outer(sign.vector, sign.vector)

    # Invert
    G_A <- tryCatch({
        solve(R_S_signed)
    }, error = function(e) {
        return(NULL)
    })

    if (is.null(G_A)) return(NULL)

    # A_A = (1' G 1)^-0.5
    val <- sum(G_A) # equivalent to t(1) %*% G %*% 1

    if (val <= 0) return(NULL)
    A_A <- val^(-0.5)

    # w_A = A_A * G * 1
    w_A <- A_A * rowSums(G_A)

    # 3. Compute 'a' vector for all variables
    w_signed <- w_A * sign.vector
    a_vec <- Rx[, active.set, drop=FALSE] %*% w_signed
    a_vec <- as.numeric(a_vec)

    # 4. Find Minimum Step Size (Gamma)
    r_j <- current.correlations[available.vars]
    a_j <- a_vec[available.vars]

    # Calculate possible gammas
    eps <- 1e-10

    # Positive sign entry
    num_plus <- r_A - r_j
    den_plus <- A_A - a_j
    g_plus <- ifelse(den_plus > eps, num_plus / den_plus, Inf)
    g_plus[g_plus < 0] <- Inf

    # Negative sign entry
    num_minus <- r_A + r_j
    den_minus <- A_A + a_j
    g_minus <- ifelse(den_minus > eps, num_minus / den_minus, Inf)
    g_minus[g_minus < 0] <- Inf

    # Find min for each var
    g_min <- pmin(g_plus, g_minus)
    s_min <- ifelse(g_plus < g_minus, 1, -1)

    # Find global winner
    best_idx_local <- which.min(g_min)

    if (length(best_idx_local) == 0 || is.infinite(g_min[best_idx_local])) {
        return(NULL)
    }

    winner_global_idx <- available.vars[best_idx_local]
    winner_gamma <- g_min[best_idx_local]
    winner_sign <- s_min[best_idx_local]

    return(list(
        next_var = winner_global_idx,
        next_sign = winner_sign,
        gamma = winner_gamma,
        a_vec = a_vec
    ))
}
