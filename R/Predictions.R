#'
#' @title Predictions for srlars Object
#'
#' @description \code{predict.srlars} returns the predictions for a srlars object.
#'
#' @method predict srlars
#'
#' @param object An object of class srlars
#' @param newx New data for predictions.
#' @param group_index Groups included in the ensemble. Default setting includes all the groups.
#' @param dynamic Argument to determine whether dynamic predictions are used based on deviating cells. Default is FALSE.
#' @param ... Additional arguments for compatibility.
#'
#' @return The predictions for the srlars object.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @seealso \code{\link{srlars}}
#'
#' @examples
#' # Required library
#' library(mvnfast)
#'
#' # Simulation parameters
#' n <- 50
#' p <- 500
#' rho.within <- 0.8
#' rho.between <- 0.2
#' p.active <- 100
#' group.size <- 25
#' snr <- 3
#' contamination.prop <- 0.2
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
#' # Contamination of data
#' contamination_indices <- 1:floor(n*contamination.prop)
#' k_lev <- 2
#' k_slo <- 100
#' x_train <- x
#' y_train <- y
#' beta_cont <- true.beta
#' beta_cont[true.beta!=0] <- beta_cont[true.beta!=0]*(1 + k_slo)
#' beta_cont[true.beta==0] <- k_slo*max(abs(true.beta))
#' for(cont_id in contamination_indices){
#'
#'   a <- runif(p, min = -1, max = 1)
#'   a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
#'   x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) +
#'     k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
#'   y_train[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
#' }
#'
#' # Ensemble models
#' ensemble_fit <- srlars(x_train, y_train,
#'                        n_models = 5,
#'                        model_saturation = c("fixed", "p-value")[1],
#'                        alpha = 0.05, model_size = n - 1,
#'                        robust = TRUE,
#'                        compute_coef = TRUE,
#'                        en_alpha = 1/4)
#'
#' # Ensemble coefficients
#' ensemble_coefs <- coef(ensemble_fit, group_index = 1:ensemble_fit$n_models)
#' sens_ensemble <- sum(which((ensemble_coefs[-1]!=0)) <= p.active)/p.active
#' spec_ensemble <- sum(which((ensemble_coefs[-1]!=0)) <= p.active)/sum(ensemble_coefs[-1]!=0)
#'
#' # Simulation of test data
#' m <- 2e3
#' x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
#' y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)
#'
#' # Prediction of test samples
#' ensemble_preds <- predict(ensemble_fit, newx = x_test,
#'                           group_index = 1:ensemble_fit$n_models,
#'                           dynamic = FALSE)
#' mspe_ensemble <- mean((y_test - ensemble_preds)^2)/sigma^2
#'
predict.srlars <- function(object, newx, group_index = NULL,
                           dynamic = FALSE,
                           ...){

  if(!dynamic){

    ensemble.coef <- coef(object, group_index = group_index)
    output <- ensemble.coef[1] + as.numeric(newx %*% ensemble.coef[-1])
    return(output)

  } else{

    DDC_cells <- cellWise::DDCpredict(newx, object$DDCx)$indcells
    x_test_cells <- newx
    x_test_cells[DDC_cells] <- NA
    cells_id <- apply(x_test_cells, 1, function(x) return(which(is.na(x))))
    var_selections <- lapply(object$coefficients, function(x) return(which(x!=0)))
    selected_models <- matrix(0, nrow = nrow(newx), ncol = object$n_models)
    for(model_id in 1:object$n_models){

      selected_models[, model_id] <- sapply(cells_id, function(x) return(any(x %in% var_selections[[model_id]])),
                                            simplify = TRUE)
    }
    selected_models <- lapply(1:nrow(newx), function(x, selected_models) return(which(selected_models[x, ] != 0)),
                              selected_models = selected_models)
    output <- numeric(nrow(newx))
    for(obs_id in 1:nrow(newx)){

      ensemble.coef <- coef(object, group_index = selected_models[[obs_id]])
      output[obs_id] <- ensemble.coef[1] + as.numeric(newx[obs_id,] %*% ensemble.coef[-1])
    }
    return(output)
  }
}
