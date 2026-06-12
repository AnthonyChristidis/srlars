#' @title Compute CV Error (Internal)
#'
#' @description
#' Evaluates the cross-validation error of a given active set using a hyper-fast 
#' Cholesky solver (with optional Huber IRLS robust fitting) and robust validation loss.
#'
#' @param cv_data List of fold data (x_train, y_train, x_val, y_val).
#' @param active_set Integer vector of active predictors.
#' @param cv_fit Character. Fitting method: "ls" or "huber" (IRLS).
#' @param cv_loss Character. Loss function: "huber", "trimmed", or "mse".
#'
#' @return Numeric. The averaged cross-validation error.
#'
#' @keywords internal
#'
#' @importFrom stats lm.fit mad quantile sd
computeCVError <- function(cv_data, active_set, cv_fit, cv_loss) {

    cv_folds <- length(cv_data)
    fold_errors <- numeric(cv_folds)

    for (f in 1:cv_folds) {
        dat <- cv_data[[f]]

        x_train <- dat$x_train
        y_train <- dat$y_train
        x_val <- dat$x_val
        y_val <- dat$y_val

        # _____________________________
        # 1. Base Case: Intercept Only
        # _____________________________
        
        if (length(active_set) == 0) {
            if (cv_fit == "huber") {
                intercept <- median(y_train)
            } else {
                intercept <- mean(y_train)
            }
            r_train <- y_train - intercept
            y_pred <- rep(intercept, length(y_val))
        } else {
            # _______________________
            # 2. Setup Data Matrices
            # _______________________

            vars <- active_set
            if (length(vars) >= nrow(x_train)) {
                vars <- vars[1:(nrow(x_train) - 1)]
            }

            X_tr <- cbind(1, x_train[, vars, drop = FALSE])
            X_te <- cbind(1, x_val[, vars, drop = FALSE])

            # ___________________________________
            # 3. Model Fitting: LS or Huber IRLS
            # ___________________________________

            if (cv_fit == "ls") {
                # --- Fast LS (Normal Equations + Ridge) ---
                XtX <- crossprod(X_tr)
                Xty <- crossprod(X_tr, y_train)
                diag(XtX) <- diag(XtX) + 1e-7

                coefs <- tryCatch({
                    solve(XtX, Xty)
                }, error = function(e) {
                    return(lm.fit(X_tr, y_train)$coefficients)
                })

            } else if (cv_fit == "huber") {
                # --- Fast 2-step IRLS for Huber M-estimator ---
                # Step A: Initial LS fit
                XtX <- crossprod(X_tr)
                diag(XtX) <- diag(XtX) + 1e-7
                coefs <- tryCatch({ solve(XtX, crossprod(X_tr, y_train)) }, 
                                  error = function(e) lm.fit(X_tr, y_train)$coefficients)
                coefs[is.na(coefs)] <- 0
                
                # Step B: 2 iterations of Reweighting
                for (iter in 1:2) {
                    resid <- y_train - as.numeric(X_tr %*% coefs)
                    scale_est <- mad(resid)
                    if (scale_est < 1e-6) scale_est <- max(sd(resid), 1e-6)
                    
                    # Huber weights: w = 1 if |resid/s| <= 1.345, else 1.345 / |resid/s|
                    u <- abs(resid / scale_est)
                    w <- ifelse(u <= 1.345, 1, 1.345 / u)
                    
                    # Weighted LS: solve (X' W X) b = X' W y
                    Xw <- X_tr * w
                    XtX_w <- crossprod(X_tr, Xw)
                    diag(XtX_w) <- diag(XtX_w) + 1e-7
                    coefs <- tryCatch({ solve(XtX_w, crossprod(Xw, y_train)) }, 
                                      error = function(e) coefs) # keep old coefs if singular
                }
            }
            
            coefs[is.na(coefs)] <- 0
            y_pred <- as.numeric(X_te %*% coefs)
            
            # Needed for Huber CV loss later
            if (cv_loss == "huber") {
                r_train <- y_train - as.numeric(X_tr %*% coefs)
            }
        }

        # _________________________________________
        # 4. Validation Residuals & Robust Scoring
        # _________________________________________

        r_val <- y_val - y_pred

        if (cv_loss == "mse") {
            fold_errors[f] <- mean(r_val^2)

        } else if (cv_loss == "trimmed") {
            r2 <- r_val^2
            cutoff <- quantile(r2, 0.90, names = FALSE)
            fold_errors[f] <- mean(r2[r2 <= cutoff])

        } else if (cv_loss == "huber") {
            k <- 1.345
            scale_est <- mad(r_train)
            if (scale_est < 1e-6) scale_est <- max(sd(r_train), 1e-6)

            u <- r_val / scale_est
            huber_loss <- ifelse(abs(u) <= k,
                                 0.5 * u^2,
                                 k * abs(u) - 0.5 * k^2)
            fold_errors[f] <- mean(huber_loss) * (scale_est^2)
        }
    }

    return(mean(fold_errors))
}