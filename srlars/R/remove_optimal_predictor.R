# -------------------------
# Remove Optimal Predictor
# -------------------------

remove_optimal_predictor <- function(lars_model, optimal_predictor,
                                     x_imp, y_imp,
                                     model_saturation, alpha){

  # Current optimal predictor
  optimal_current <- lars_model[["optimal_predictor"]]

  # Remove optimal predictor
  position_optimal <- c(which(lars_model[["candidate_predictors"]] == optimal_predictor))
  lars_model[["candidate_predictors"]] <-
    lars_model[["candidate_predictors"]][-position_optimal]

  # Adjust a_vec
  lars_model[["a_vec"]] <-
    lars_model[["a_vec"]][-position_optimal]

  # Adjust r_vec
  lars_model[["r_vec"]] <-
    lars_model[["r_vec"]][-position_optimal]

  # Adjust gamma_vec_p
  lars_model[["gamma_vec_p"]] <-
    lars_model[["gamma_vec_p"]][-position_optimal]

  # Adjust gamma_vec_m
  lars_model[["gamma_vec_m"]] <-
    lars_model[["gamma_vec_m"]][-position_optimal]

  # Adjust gamma_vec
  lars_model[["gamma_vec"]] <-
    apply(cbind(lars_model[["gamma_vec_p"]], lars_model[["gamma_vec_m"]]), 1, min)

  # Adjust gamma
  lars_model[["gamma"]] <-
    min(lars_model[["gamma_vec"]])

  # Adjust optimal predictor
  lars_model[["optimal_predictor"]] <-
    lars_model[["candidate_predictors"]][which.min(lars_model[["gamma_vec"]])]

  # Adjustment of model fit if new optimal predictor
  if(lars_model[["optimal_predictor"]] != optimal_current){

    # Update candidate MM model
    lars_model[["rss_candidate"]] <-
      mean((y_imp - x_imp[, c(lars_model[["model_predictors"]], lars_model[["optimal_predictor"]])] %*%
              solve(t(x_imp[, c(lars_model[["model_predictors"]], lars_model[["optimal_predictor"]])]) %*%
                      x_imp[, c(lars_model[["model_predictors"]], lars_model[["optimal_predictor"]])]) %*%
              t(x_imp[, c(lars_model[["model_predictors"]], lars_model[["optimal_predictor"]])]) %*% y_imp)^2)

    # Update p-value
    lars_model[["p_value"]] <-
      pf((lars_model[["rss_current"]] - lars_model[["rss_candidate"]]) *
           (nrow(x_imp) - length(lars_model[["model_predictors"]]) - 1) /
           lars_model[["rss_candidate"]], df1 = 1, df2 = nrow(x_imp) - length(lars_model[["model_predictors"]] - 1),
         lower.tail = TRUE)

    # Update model saturation
    if((model_saturation == "p-value") && (lars_model[["p_value"]] >= alpha))
      lars_model[["saturated"]] <- TRUE
  }

  return(lars_model)
}
