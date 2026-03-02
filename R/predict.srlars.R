#' @title Predictions for srlars Object
#'
#' @description \code{predict.srlars} returns the predictions for a srlars object.
#'
#' @method predict srlars
#'
#' @param object An object of class srlars.
#' @param newx New data matrix for predictions.
#' @param model_index Indices of the sub-models to include in the ensemble. Default is NULL (all models).
#' @param dynamic Logical. If TRUE, and the model was trained robustly, the new data \code{newx} is cleaned using
#'        \code{\link[cellWise]{DDCpredict}} before prediction. This ensures consistency with the robust training phase.
#'        Default is TRUE.
#' @param ... Additional arguments for compatibility.
#'
#' @return A numeric vector of predictions.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @seealso \code{\link{srlars}}
#'
#' @importFrom cellWise DDCpredict
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
#'                        tolerance = 0.01,
#'                        robust = TRUE,
#'                        compute_coef = TRUE)
#'
#' # Simulation of test data
#' m <- 50
#' x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
#' y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)
#'
#' # Prediction of test samples
#' # Default: Averaged prediction from all models
#' # Dynamic imputation is used by default if the model is robust
#' ensemble_preds <- predict(ensemble_fit, newx = x_test)
#'
#' # Evaluate MSPE
#' mspe <- mean((y_test - ensemble_preds)^2) / sigma^2
#' print(paste("MSPE:", mspe))
#'
predict.srlars <- function(object,
                           newx,
                           model_index = NULL,
                           dynamic = TRUE,
                           ...) {

    # 1. Validate Inputs
    newx <- as.matrix(newx)

    if (is.null(colnames(newx))) {
        colnames(newx) <- paste0("V", 1:ncol(newx))
    }

    if(is.null(model_index)){
        model_index <- 1:object$n_models
    } else{
        if(any(!(model_index %in% 1:object$n_models)))
            stop("The model_index contains invalid indices.")
    }

    # 2. Dynamic Imputation (Robust Prediction)
    x_for_pred <- newx

    # Check if we should (and can) perform robust cleaning
    if (dynamic && isTRUE(object$robust) && !is.null(object$ddc.object)) {

        # We use DDCpredict to clean the new data.
        # DDCpredict requires the input to have the same number of columns as the training data.
        # Since training data was [X, y] (p+1 columns), we must augment newx with a dummy response column.

        # Append dummy column of NAs
        newx_aug <- cbind(newx, NA)

        # Ensure column names match expected pattern if possible.
        # We append a placeholder name for the dummy response.
        colnames(newx_aug) <- c(colnames(newx), "response_placeholder")

        ddc_pred <- tryCatch({
            # Run DDCpredict on augmented data
            # Note: We rely on DDCpredict to ignore the NA column for cleaning X, or treat it as missing.
            cellWise::DDCpredict(Xnew = newx_aug, InitialDDC = object$ddc.object)
        }, error = function(e) {
            warning(paste("DDCpredict failed:", e$message, "Falling back to raw newx."))
            return(NULL)
        })

        if (!is.null(ddc_pred)) {
            # Extract the cleaned X part (remove the dummy response column)
            # DDCpredict returns $Ximp with same dimensions as input
            x_for_pred <- ddc_pred$Ximp[, 1:ncol(newx), drop = FALSE]
        }
    }

    # 3. Compute Predictions
    n_test <- nrow(x_for_pred)
    final_preds <- numeric(n_test)
    n_groups <- length(model_index)

    for (k in model_index) {
        # Extract params
        beta_k <- object$coefficients[[k]]
        intercept_k <- object$intercepts[k]

        # Linear predictor: alpha + X * beta
        preds_k <- intercept_k + (x_for_pred %*% beta_k)

        # Accumulate
        final_preds <- final_preds + preds_k
    }

    # Average
    if (n_groups > 0) {
        final_preds <- final_preds / n_groups
    }

    return(as.numeric(final_preds))
}
