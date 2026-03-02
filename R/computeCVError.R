#' @title Compute Cross-Validated Prediction Error (Internal)
#'
#' @description
#' Computes the V-fold CV error for a linear regression model using fast LS fitting.
#'
#' @param x Design matrix (imputed).
#' @param y Response vector (imputed).
#' @param active_set Indices of predictors to include.
#' @param folds Number of folds (default 5).
#'
#' @return A single numeric value (Mean Squared Prediction Error).
#'
#' @keywords internal
#'
#' @importFrom stats lm.fit
#'
computeCVError <- function(x, y,
                           active_set,
                           folds = 5) {

    n <- nrow(x)

    # 1. Base Case: Intercept Only
    if (length(active_set) == 0) {
        # Prediction is mean(y). CV error approx variance.
        # More precise: sum((y - mean(y))^2) / n
        return(mean((y - mean(y))^2))
    }

    # 2. Setup Data
    x_sub <- x[, active_set, drop = FALSE]
    # Add intercept manually.
    x_design <- cbind(1, x_sub)

    fold_ids <- sample(rep(1:folds, length.out = n))
    total_sq_error <- 0

    # 3. Fast CV Loop
    for (k in 1:folds) {
        test_idx <- which(fold_ids == k)
        train_idx <- which(fold_ids != k)

        # Train Data
        X_tr <- x_design[train_idx, , drop=FALSE]
        y_tr <- y[train_idx]

        # Test Data
        X_te <- x_design[test_idx, , drop=FALSE]

        # --- OPTIMIZED SOLVER (Normal Equations) ---
        # Solve: (X'X) * beta = X'y

        # 1. Compute Cross Products (very fast in R)
        XtX <- crossprod(X_tr)
        Xty <- crossprod(X_tr, y_tr)

        # 2. Add Ridge for Stability (The "Safety Net")
        # If matrix is singular, this prevents crash.
        # 1e-7 is small enough to not bias result but big enough for Cholesky.
        diag(XtX) <- diag(XtX) + 1e-7

        # 3. Solve using Cholesky
        coefs <- tryCatch({
            solve(XtX, Xty)
        }, error = function(e) {
            # Fallback to QR (lm.fit) only if Cholesky fails completely
            return(lm.fit(X_tr, y_tr)$coefficients)
        })

        # Predict
        y_pred <- X_te %*% coefs

        # Accumulate
        total_sq_error <- total_sq_error + sum((y[test_idx] - y_pred)^2)
    }

    return(total_sq_error / n)
}
