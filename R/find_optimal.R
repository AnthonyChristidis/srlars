# -----------------------------------
# Function to Find Optimal Predictor
# -----------------------------------

find_optimal <- function(lars_model,
                         x_imp, y_imp,
                         Rx, Ry,
                         model_saturation,
                         alpha){

  # Model saturation check
  if(det(Rx[lars_model[["model_predictors"]], lars_model[["model_predictors"]], drop = FALSE]) < 1e-8){

    lars_model[["saturated"]] <- TRUE
    return(lars_model)
  }

  # Find inverse of matrix needed for computation
  DRD_inv <-
    solve(diag(lars_model[["s_vec"]], nrow = length(lars_model[["s_vec"]])) %*%
            Rx[lars_model[["model_predictors"]], lars_model[["model_predictors"]]] %*%
            diag(lars_model[["s_vec"]], nrow = length(lars_model[["s_vec"]])))

  # Update a
  lars_model[["a"]] <-
    as.numeric((t(rep(1, length(lars_model[["model_predictors"]]))) %*%
                  DRD_inv %*%
                  rep(1, length(lars_model[["model_predictors"]])))^(-1/2))

  # Update w_vec
  lars_model[["w_vec"]] <-
    lars_model[["a"]] * DRD_inv %*% rep(1, length(lars_model[["model_predictors"]]))

  # Update a_vec
  lars_model[["a_vec"]] <-
    c(t(diag(lars_model[["s_vec"]], nrow = length(lars_model[["s_vec"]])) %*%
          Rx[lars_model[["model_predictors"]], lars_model[["candidate_predictors"]], drop = FALSE]) %*%
        lars_model[["w_vec"]])

  # Update gamma_vec_p
  lars_model[["gamma_vec_p"]] <-
    (lars_model[["r"]] - lars_model[["r_vec"]]) /
    (lars_model[["a"]] - lars_model[["a_vec"]])

  # Update gamma_vec_m
  lars_model[["gamma_vec_m"]] <-
    (lars_model[["r"]] + lars_model[["r_vec"]]) /
    (lars_model[["a"]] + lars_model[["a_vec"]])

  # Update gamma_vec
  lars_model[["gamma_vec"]] <-
    apply(cbind(lars_model[["gamma_vec_p"]], lars_model[["gamma_vec_m"]]), 1, min)

  # Update gamma
  lars_model[["gamma"]] <-
    min(lars_model[["gamma_vec"]])

  # Update optimal predictor
  lars_model[["optimal_predictor"]] <-
    lars_model[["candidate_predictors"]][which.min(lars_model[["gamma_vec"]])]

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

  return(lars_model)
}
