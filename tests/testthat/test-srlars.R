test_that("srlars runs correctly on synthetic data", {

  # --- Setup: Generate Simple Data ---
  set.seed(0)
  n <- 50
  p <- 20
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("V", 1:p) # Ensure names exist
  beta <- c(rep(2, 3), rep(0, p - 3))
  y <- as.numeric(x %*% beta + rnorm(n))

  # --- Test 1: Robust Execution (DDC + wrap correlations) ---
  fit <- srlars(
    x, y,
    n_models = 3,
    tolerance = 1e-8,
    x_preprocess = "ddc",
    y_preprocess = "wrap",
    cor_estimator = "wrap",
    cv_preprocess = "global",
    cv_fit = "huber",
    cv_loss = "huber",
    cv_folds = 5,
    compute_coef = TRUE
  )

  # Check Object Structure
  expect_s3_class(fit, "srlars")
  expect_equal(fit$n_models, 3)
  expect_equal(fit$x_preprocess, "ddc")
  expect_equal(fit$y_preprocess, "wrap")
  expect_type(fit$active.sets, "list")
  expect_type(fit$coefficients, "list")

  # Check DDC object presence (needed for predict if implemented)
  expect_false(is.null(fit$ddc.object))

  # --- Test 2: Non-Robust Execution ---
  fit_nr <- srlars(
    x, y,
    n_models = 2,
    x_preprocess = "none",
    y_preprocess = "none",
    cor_estimator = "pearson",
    cv_preprocess = "global",
    cv_fit = "ls",
    cv_loss = "mse",
    cv_folds = 5,
    compute_coef = TRUE
  )

  expect_s3_class(fit_nr, "srlars")
  expect_equal(fit_nr$x_preprocess, "none")
  expect_equal(fit_nr$y_preprocess, "none")
  expect_true(is.null(fit_nr$ddc.object))

  # --- Test 3: Coefficients Method ---
  # Should return vector of length p + 1 (intercept)
  coefs <- coef(fit)
  expect_type(coefs, "double")
  expect_length(coefs, p + 1)

  # --- Test 4: Prediction Method ---
  newx <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
  colnames(newx) <- colnames(x)

  # Robust prediction
  # (Only pass dynamic=TRUE if your predict.srlars supports it.)
  preds <- predict(fit, newx)
  expect_type(preds, "double")
  expect_length(preds, 10)

  # Non-robust prediction
  preds_nr <- predict(fit_nr, newx)
  expect_type(preds_nr, "double")
  expect_length(preds_nr, 10)

  # --- Test 5: Input Checks ---
  expect_error(srlars(x, y[1:(n - 1)]))
})