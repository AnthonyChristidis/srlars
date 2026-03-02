#' @title Coefficients for srlars Object
#'
#' @description
#' \code{coef.srlars} returns the averaged coefficients for a srlars object.
#'
#' @method coef srlars
#'
#' @param object An object of class srlars.
#' @param model_index Indices of the sub-models to include in the ensemble average. Default is NULL, which includes all models.
#' @param ... Additional arguments for compatibility.
#'
#' @return A numeric vector containing the averaged intercept (first element) and slope coefficients.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @seealso \code{\link{srlars}}, \code{\link{predict.srlars}}
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
#' # Ensemble coefficients
#' # Default: Average over all models
#' ensemble_coefs <- coef(ensemble_fit)
#'
#' # Sensitivity (Recall)
#' active_selected <- which(ensemble_coefs[-1] != 0)
#' true_active <- which(true.beta != 0)
#' recall <- length(intersect(active_selected, true_active)) / length(true_active)
#' print(paste("Recall:", recall))
#'
#' # Precision
#' if(length(active_selected) > 0){
#'   precision <- length(intersect(active_selected, true_active)) / length(active_selected)
#' } else {
#'   precision <- 0
#' }
#' print(paste("Precision:", precision))
#'
coef.srlars <- function(object, model_index = NULL, ...) {

    # 1. Validate Group Index
    if (is.null(model_index)) {
        model_index <- 1:object$n_models
    } else {
        if (any(model_index < 1) || any(model_index > object$n_models)) {
            stop("The model_index contains invalid indices.")
        }
    }

    # 2. Determine p
    # The first coefficient vector gives us the dimension
    p <- length(object$coefficients[[1]])

    # 3. Initialize Accumulator
    # Length p + 1 for Intercept + Slopes
    accumulated_coefs <- numeric(p + 1)

    # 4. Average Coefficients
    n_groups <- length(model_index)

    for (k in model_index) {
        # Extract current model's params
        beta_k <- object$coefficients[[k]]
        intercept_k <- object$intercepts[k]

        # Add to accumulator (Intercept, then Betas)
        accumulated_coefs <- accumulated_coefs + c(intercept_k, beta_k)
    }

    # Divide by number of groups to get the mean
    final_coef <- accumulated_coefs / n_groups

    # 5. Add Names (Optional but nice)
    # names(final_coef) <- c("(Intercept)", paste0("V", 1:p))

    return(final_coef)
}

