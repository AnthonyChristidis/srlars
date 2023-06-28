# ----------------------
# Add Optimal Predictor
# ----------------------

add_optimal_predictor <- function(lars_model,
                                  model_saturation, model_size, n){

  # Add optimal predictor
  lars_model[["model_predictors"]] <-
    c(lars_model[["model_predictors"]], lars_model[["optimal_predictor"]])

  # Remove optimal predictor from candidates
  position_optimal <- c(which(lars_model[["candidate_predictors"]] == lars_model[["optimal_predictor"]]))
  lars_model[["candidate_predictors"]] <-
    lars_model[["candidate_predictors"]][-position_optimal]

  # Update s_vec
  lars_model[["s_vec"]] <-
    c(lars_model[["s_vec"]],
      ifelse(lars_model[["gamma_vec_p"]][position_optimal] <
               lars_model[["gamma_vec_m"]][position_optimal], 1, -1))

  # Update r
  lars_model[["r"]] <-
    lars_model[["r"]] - lars_model[["gamma"]] * lars_model[["a"]]

  # Update r_vec
  lars_model[["r_vec"]] <-
    lars_model[["r_vec"]][-position_optimal] - lars_model[["gamma"]] * lars_model[["a_vec"]][-position_optimal]

  # Check if model saturation is achieved
  if(length(lars_model[["model_predictors"]]) == (n - 1))
    lars_model[["saturated"]] <- TRUE else if(model_saturation == "fixed" && length(lars_model[["model_predictors"]]) == (n - 1))
      lars_model[["saturated"]] <- TRUE

  # Update current MM model
  lars_model[["rss_current"]] <-lars_model[["rss_candidate"]]

  return(lars_model)
}
