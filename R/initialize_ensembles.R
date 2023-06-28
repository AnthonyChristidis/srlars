# ---------------------------------
# Function to Initialize Ensembles
# ---------------------------------

initialize_ensembles <- function(lars_models,
                                 x_imp, y_imp,
                                 n_models, Ry){

  initial_predictors <- order(abs(Ry), decreasing = TRUE)[1:n_models]

  for(model_id in 1:n_models){

    lars_models[[model_id]] <- list()
    lars_models[[model_id]][["model_predictors"]] <- initial_predictors[model_id]
    lars_models[[model_id]][["candidate_predictors"]] <- c(1:length(Ry))[-c(initial_predictors)]
    lars_models[[model_id]][["s_vec"]] <- sign(Ry[initial_predictors[model_id]])
    lars_models[[model_id]][["a"]] <- numeric(0)
    lars_models[[model_id]][["w_vec"]] <- numeric(0)
    lars_models[[model_id]][["a_vec"]] <- numeric(0)
    lars_models[[model_id]][["r"]] <- Ry[initial_predictors[model_id]]
    lars_models[[model_id]][["r_vec"]] <- Ry[-c(initial_predictors)]
    lars_models[[model_id]][["gamma"]] <- numeric(0)
    lars_models[[model_id]][["gamma_vec"]] <- numeric(0)
    lars_models[[model_id]][["gamma_vec_p"]] <- numeric(0)
    lars_models[[model_id]][["gamma_vec_m"]] <- numeric(0)
    lars_models[[model_id]][["optimal_predictor"]] <- numeric(0)
    lars_models[[model_id]][["rss_current"]] <-
      mean((y_imp - x_imp[, lars_models[[model_id]][["model_predictors"]]] %*%
              solve(t(x_imp[, lars_models[[model_id]][["model_predictors"]]]) %*%
                      x_imp[, lars_models[[model_id]][["model_predictors"]]]) %*%
              t(x_imp[, lars_models[[model_id]][["model_predictors"]]]) %*% y_imp)^2)
    lars_models[[model_id]][["rss_candidate"]] <- numeric(0)
    lars_models[[model_id]][["p_value"]] <- numeric(0)
    lars_models[[model_id]][["saturated"]] <- FALSE
  }

  return(lars_models)
}
