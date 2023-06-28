#'
#' @importFrom stats coef sd median mad cor pf
#'
#' @title Robust Split Least Angle Regression
#'
#' @description \code{srlars} performs split robust least angle regression.
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param n_models Number of models into which the variables are split.
#' @param model_saturation Criterion to determine if a model is saturated. Must be one of "fixed" (default) or "p-value".
#' @param alpha P-value used to determine when the model is saturated
#' @param model_size Size of the models in the ensemble.
#' @param robust Argument to determine if robust measures of location, scale and correlation are used. Default is TRUE.
#' @param compute_coef Argument to determine if coefficients are computed (via the elastic net) for each model. Default is FALSE.
#' @param en_alpha Elastic net mixing parmeter for parameters shrinkage. Default is 1/4.
#'
#' @return An object of class srlars
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @seealso \code{\link{coef.srlars}}, \code{\link{predict.srlars}}
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
srlars <- function(x, y,
                   n_models = 1,
                   model_saturation = c("fixed", "p-value")[1],
                   alpha = 0.05,
                   model_size = NULL,
                   robust = TRUE,
                   compute_coef = FALSE,
                   en_alpha = 1/4){


  # Data input check
  DataCheck(x, y,
            n_models,
            model_saturation,
            alpha,
            model_size,
            robust,
            compute_coef,
            en_alpha)

  # Shuffle the data
  n <- nrow(x)
  p <- ncol(x)
  random.permutation <- sample(1:n, n)
  x <- x[random.permutation, ]
  y <- y[random.permutation]

  # Model size check
  if(model_saturation == "fixed"){

    if(is.null(model_size))
      model_size <- n - 1 else if(model_size >= n)
        stop("The size of the models cannot be equal or exceed n.")
  }

  # Standarization predictors and response
  # Computation of correlation for predictors and response
  if(robust){

    # Standardization of predictors and response
    x.std <- apply(x, 2, function(x) return((x - median(x))/mad(x)))
    y.std <- (y - median(y))/mad(y)

    # Computation of correlations for predictors and response
    xy.std <- cbind(x.std, y.std)
    DDCxy <- cellWise::DDC(xy.std, DDCpars = list(fastDDC = TRUE, silent = TRUE))
    rob.cor <- cor(DDCxy$Ximp)
    Rx <- rob.cor[-nrow(rob.cor),-ncol(rob.cor)]
    Ry <- rob.cor[-nrow(rob.cor), ncol(rob.cor)]

    # Full data
    DDCxy <- cellWise::DDC(cbind(x, y), DDCpars = list(fastDDC = TRUE, silent = TRUE))
    x_imp <- DDCxy$Ximp[, -ncol(DDCxy$Ximp)]
    y_imp <- DDCxy$Ximp[, ncol(DDCxy$Ximp)]

  } else{

    # Standardization of predictors and response
    x.std <- apply(x, 2, function(x) return((x - mean(x))/sd(x)))
    y.std <- (y - mean(y))/sd(y)

    # Computation of correlations for predictors and response
    xy.std <- cbind(x.std, y.std)
    CORxy <- cor(xy.std)
    Rx <- CORxy[-nrow(CORxy), -ncol(CORxy)]
    Ry <- CORxy[-nrow(CORxy), ncol(CORxy)]

    # Full data
    x_imp <- x
    y_imp <- y
  }

  # Initialize LARS models
  lars_models <- list()
  lars_models <- initialize_ensembles(lars_models,
                                      x_imp, y_imp,
                                      n_models, Ry)

  # Vector of unsaturated models
  unsaturated_models <- 1:n_models

  # Vector of p-values
  p_values <- numeric(5)

  # Finding optimal predictor
  for(model_id in unsaturated_models)
    lars_models[[model_id]] <- find_optimal(lars_models[[model_id]],
                                            x_imp, y_imp,
                                            Rx, Ry,
                                            model_saturation, alpha)

  # Looping over the models
  while(length(unsaturated_models) != 0){

    # Finding optimal unsaturated model
    p_values <- sapply(lars_models, function(x) x$p_value)
    p_values[-c(unsaturated_models)] <- Inf
    optimal_model <- which.min(p_values)
    optimal_predictor <- lars_models[[optimal_model]][["optimal_predictor"]]

    # Add optimal predictor for optimal model
    lars_models[[optimal_model]] <- add_optimal_predictor(lars_models[[optimal_model]],
                                                          model_saturation, model_size, n)
    # Update optimal model
    if(!lars_models[[optimal_model]][["saturated"]])
      lars_models[[optimal_model]] <- find_optimal(lars_models[[optimal_model]],
                                                   x_imp, y_imp,
                                                   Rx, Ry,
                                                   model_saturation, alpha)

    # Remove optimal predictor for non-optimal unsaturated models
    for(model_id in unsaturated_models[-c(which(unsaturated_models == optimal_model))])
      lars_models[[model_id]] <- remove_optimal_predictor(lars_models[[model_id]], optimal_predictor,
                                                          x_imp, y_imp,
                                                          model_saturation, alpha)

    # Update unsaturated models
    unsaturated_models <- which(!sapply(lars_models, function(x) x$saturated))
  }

  # Selected predictors
  selections <- lapply(lars_models, function(x) x$model_predictors)

  # Creating the list for the output
  output <- list(x = x, y = y,
                 n_models = n_models,
                 model_saturation = model_saturation,
                 alpha = alpha,
                 model_size = model_size,
                 robust = robust,
                 compute_coef = compute_coef,
                 selections = selections,
                 intercepts = list(), coefficients = list(),
                 DDCx = cellWise::DDC(x, DDCpars = list(fastDDC = TRUE, silent = TRUE)))

  # Computation of final coefficients
  if(compute_coef){

    for(model_id in 1:n_models){

      en_fit <- glmnet::cv.glmnet(x_imp[, output$selections[[model_id]]], y_imp,
                                  alpha = en_alpha)
      output$intercepts[[model_id]] <- coef(en_fit, s = "lambda.min")[1]
      output$coefficients[[model_id]] <- numeric(p)
      output$coefficients[[model_id]][output$selections[[model_id]]] <- coef(en_fit, s = "lambda.min")[-1]
    }
  }

  # Create the object of class "srlars"
  class(output) <- append("srlars", class(output))

  # Returning the output from the stepwise algorithm
  return(output)
}



