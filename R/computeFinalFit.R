#' @title FSCRE Final Robust Fitting (Internal)
#'
#' @description
#' Implements Stage 3: Fitting robust MM-estimators to the selected variable sets.
#'
#' @param x.imp Imputed design matrix.
#' @param y.imp Imputed response vector.
#' @param active.sets List of integer vectors (indices of selected variables).
#' @param compute_coef Logical.
#'
#' @return A list containing `coefficients` (list of vectors) and `intercepts` (vector).
#'
#' @keywords internal
#'
#' @importFrom robustbase lmrob
#' @importFrom stats lm.fit coef median
#'
computeFinalFit <- function(x.imp, y.imp,
                            active.sets,
                            compute_coef) {

    n_models <- length(active.sets)
    p <- ncol(x.imp)
    n <- nrow(x.imp)

    # Initialize outputs
    intercepts <- numeric(n_models)
    coefficients <- vector("list", n_models)

    # If not computing coefficients, return empty structures
    if (!compute_coef) {
        for(k in 1:n_models) coefficients[[k]] <- numeric(p)
        return(list(coefficients = coefficients, intercepts = intercepts))
    }

    # Robust fitting loop
    for (k in 1:n_models) {
        vars <- active.sets[[k]]

        # Setup default (empty) coefficient vector
        full_coefs <- numeric(p)

        if (length(vars) > 0) {

            # Constraint Check: p_k < n
            if (length(vars) >= n) {
                warning(paste("Model", k, "has", length(vars), "variables, but n =", n,
                              ". Truncating to n-1 variables for final fit."))
                vars <- vars[1:(n-1)]
            }

            # Prepare data for lmrob (requires formula interface for safety)
            dat_robust <- data.frame(y = y.imp, x.imp[, vars, drop = FALSE])

            # Fit MM-estimator with fallback
            fit <- tryCatch({
                suppressWarnings(robustbase::lmrob(y ~ ., data = dat_robust))
            }, error = function(e) {
                warning(paste("MM-fit failed for Model", k, "- falling back to LS."))

                # FALLBACK: Use lm.fit for speed
                # Manually add intercept column
                x_ls <- cbind(1, x.imp[, vars, drop = FALSE])
                return(lm.fit(x = x_ls, y = y.imp))
            })

            # Extract results based on object type
            if (inherits(fit, "lmrob")) {
                all_coefs <- coef(fit)
            } else {
                # lm.fit returns a list, not an S3 object with coef method
                all_coefs <- fit$coefficients
            }

            # First element is intercept, rest are slopes
            intercepts[k] <- all_coefs[1]
            fit_coefs <- all_coefs[-1]

            # Fill in the non-zero coefficients
            full_coefs[vars] <- fit_coefs

        } else {
            # Empty model case
            intercepts[k] <- median(y.imp)
        }

        coefficients[[k]] <- full_coefs
    }

    return(list(
        coefficients = coefficients,
        intercepts = intercepts
    ))
}
