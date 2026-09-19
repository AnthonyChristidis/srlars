test_that("cv.srlars selects max_share and dispatches coef/predict correctly", {

  set.seed(0)
  n <- 60
  p <- 15
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("V", 1:p)
  beta <- c(rep(2, 3), rep(0, p - 3))
  y <- as.numeric(x %*% beta + rnorm(n))

  cv_fit <- cv.srlars(
    x, y,
    n_models = 3,
    outer_folds = 3,
    tolerance = 1e-4,
    x_preprocess = "ddc",
    y_preprocess = "wrap",
    cor_estimator = "wrap",
    cv_preprocess = "global",
    cv_fit = "huber",
    cv_loss = "huber",
    cv_folds = 3,
    compute_coef = TRUE
  )

  # --- Structure ---
  expect_s3_class(cv_fit, "cv.srlars")
  expect_s3_class(cv_fit, "srlars")
  expect_true(cv_fit$max_share %in% cv_fit$share_grid)
  expect_true(cv_fit$max_share >= 1 && cv_fit$max_share <= 3)
  expect_length(cv_fit$cv_errors, length(cv_fit$share_grid))
  expect_equal(cv_fit$share_grid, 1:3)

  # --- Dispatch to coef.srlars / predict.srlars works without new methods ---
  coefs <- coef(cv_fit)
  expect_type(coefs, "double")
  expect_length(coefs, p + 1)

  newx <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
  colnames(newx) <- colnames(x)
  preds <- predict(cv_fit, newx)
  expect_type(preds, "double")
  expect_length(preds, 10)

  # --- Validation errors ---
  expect_error(cv.srlars(x, y, n_models = 3, share_grid = c(0, 1, 2)))
  expect_error(cv.srlars(x, y, n_models = 3, share_grid = c(1, 4)))
  expect_error(cv.srlars(x, y, n_models = 3, share_grid = c(1, 1, 2)))
  expect_error(cv.srlars(x, y, n_models = 3, outer_folds = 1))
})

test_that("computeRobustLoss extraction preserves computeCVError behavior", {

  set.seed(1)
  n <- 30
  cv_data <- list(list(
    x_train = matrix(rnorm(n * 5), n, 5),
    y_train = rnorm(n),
    x_val   = matrix(rnorm(n * 5), n, 5),
    y_val   = rnorm(n)
  ))

  for (loss in c("mse", "trimmed", "huber")) {
    for (fitm in c("ls", "huber")) {
      err <- computeCVError(cv_data, active_set = c(1, 2), cv_fit = fitm, cv_loss = loss)
      expect_true(is.numeric(err) && length(err) == 1 && is.finite(err))
    }
  }
})
