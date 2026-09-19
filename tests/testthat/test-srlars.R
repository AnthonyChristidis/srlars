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

test_that("max_share defaults to fully disjoint models and is respected when relaxed", {

  set.seed(0)
  n <- 50
  p <- 20
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("V", 1:p)
  beta <- c(rep(2, 3), rep(0, p - 3))
  y <- as.numeric(x %*% beta + rnorm(n))

  common_args <- list(
    x = x, y = y,
    n_models = 3,
    tolerance = 1e-4,
    x_preprocess = "ddc",
    y_preprocess = "wrap",
    cor_estimator = "wrap",
    cv_preprocess = "global",
    cv_fit = "huber",
    cv_loss = "huber",
    cv_folds = 5,
    compute_coef = TRUE
  )

  # --- Default equivalence: omitting max_share == max_share = 1 ---
  set.seed(123)
  fit_default <- do.call(srlars, common_args)

  set.seed(123)
  fit_explicit_1 <- do.call(srlars, c(common_args, list(max_share = 1)))

  expect_identical(fit_default$active.sets, fit_explicit_1$active.sets)
  expect_equal(fit_default$coefficients, fit_explicit_1$coefficients)
  expect_equal(fit_default$intercepts, fit_explicit_1$intercepts)

  # --- Disjointness at default ---
  sets <- fit_default$active.sets
  for (i in seq_along(sets)) {
    for (j in seq_along(sets)) {
      if (i < j) {
        expect_length(intersect(sets[[i]], sets[[j]]), 0)
      }
    }
  }

  # --- Cap respected when relaxed ---
  set.seed(123)
  fit_shared <- do.call(srlars, c(common_args, list(max_share = 3)))
  usage <- table(unlist(fit_shared$active.sets))
  expect_true(all(usage <= 3))

  # --- Validation errors ---
  expect_error(do.call(srlars, c(common_args, list(max_share = 0))))
  expect_error(do.call(srlars, c(common_args, list(max_share = 1.5))))
  expect_error(do.call(srlars, c(common_args, list(max_share = 4)))) # n_models = 3
})

test_that("max_share forces distinct seeds when 1 < max_share < n_models, but not at n_models", {

  set.seed(0)
  n <- 60
  p <- 30
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("V", 1:p)
  beta <- c(rep(2, 3), rep(0, p - 3))
  y <- as.numeric(x %*% beta + rnorm(n))

  common_args <- list(
    x = x, y = y,
    n_models = 6,
    tolerance = 1e-4,
    x_preprocess = "ddc",
    y_preprocess = "wrap",
    cor_estimator = "wrap",
    cv_preprocess = "global",
    cv_fit = "huber",
    cv_loss = "huber",
    cv_folds = 5,
    compute_coef = FALSE
  )

  first_picks <- function(active.sets) {
    vapply(active.sets, function(s) if (length(s) > 0) s[1] else NA_integer_, integer(1))
  }

  # --- 1 < max_share < n_models: seeds must be pairwise distinct ---
  set.seed(123)
  fit_mid <- do.call(srlars, c(common_args, list(max_share = 3)))
  seeds <- first_picks(fit_mid$active.sets)
  seeds_present <- seeds[!is.na(seeds)]
  expect_equal(length(seeds_present), length(unique(seeds_present)))

  # --- max_share = n_models: seed restriction lifted (no assertion on distinctness,
  #     but the pool computation must skip the seed.vars exclusion; verified by checking
  #     the run completes and respects the usage cap only) ---
  set.seed(123)
  fit_full <- do.call(srlars, c(common_args, list(max_share = 6)))
  usage_full <- table(unlist(fit_full$active.sets))
  expect_true(all(usage_full <= 6))
})

test_that("n_min defaults to no floor and forces growth of otherwise-small models when set", {

  # Pure noise: x and y are unrelated, so the normal positive-benefit/tolerance
  # stopping rule should halt selection almost immediately without n_min.
  set.seed(0)
  n <- 50
  p <- 20
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("V", 1:p)
  y <- rnorm(n)

  common_args <- list(
    x = x, y = y,
    n_models = 4,
    max_share = 1,
    tolerance = 1e-4,
    x_preprocess = "ddc",
    y_preprocess = "wrap",
    cor_estimator = "wrap",
    cv_preprocess = "global",
    cv_fit = "huber",
    cv_loss = "huber",
    cv_folds = 5,
    compute_coef = FALSE
  )

  # --- Default equivalence: omitting n_min == n_min = NULL ---
  set.seed(321)
  fit_default <- do.call(srlars, common_args)

  set.seed(321)
  fit_explicit_null <- do.call(srlars, c(common_args, list(n_min = NULL)))

  expect_identical(fit_default$active.sets, fit_explicit_null$active.sets)

  # --- Without n_min, active sets on pure noise should generally stay small ---
  sizes_default <- vapply(fit_default$active.sets, length, integer(1))

  # --- With n_min, every model must reach the floor (plenty of noise variables
  #     available, max_share = 1 disjoint, floor well within min(n - 1, p)) ---
  set.seed(321)
  fit_floor <- do.call(srlars, c(common_args, list(n_min = 3)))
  sizes_floor <- vapply(fit_floor$active.sets, length, integer(1))
  expect_true(all(sizes_floor >= 3))

  # The floor should have forced at least some growth relative to the unforced run
  expect_true(sum(sizes_floor) > sum(sizes_default))

  # --- max_share usage cap still respected even while forcing ---
  usage_floor <- table(unlist(fit_floor$active.sets))
  expect_true(all(usage_floor <= 1))

  # --- Validation errors ---
  expect_error(do.call(srlars, c(common_args, list(n_min = 0))))
  expect_error(do.call(srlars, c(common_args, list(n_min = 2.5))))
  expect_error(do.call(srlars, c(common_args, list(n_min = n)))) # exceeds n - 1
  expect_error(do.call(srlars, c(common_args, list(n_min = p + 1)))) # exceeds p
})

test_that("n_min composes with max_share seed-diversity when both are active", {

  set.seed(0)
  n <- 60
  p <- 25
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("V", 1:p)
  y <- rnorm(n) # pure noise, forces the floor to bind

  set.seed(321)
  fit <- srlars(x, y,
               n_models = 5,
               max_share = 2,
               n_min = 2,
               tolerance = 1e-4,
               x_preprocess = "ddc",
               y_preprocess = "wrap",
               cor_estimator = "wrap",
               cv_preprocess = "global",
               cv_fit = "huber",
               cv_loss = "huber",
               cv_folds = 5,
               compute_coef = FALSE)

  sizes <- vapply(fit$active.sets, length, integer(1))
  expect_true(all(sizes >= 2))

  usage <- table(unlist(fit$active.sets))
  expect_true(all(usage <= 2))

  first_picks <- vapply(fit$active.sets, function(s) if (length(s) > 0) s[1] else NA_integer_, integer(1))
  seeds_present <- first_picks[!is.na(first_picks)]
  expect_equal(length(seeds_present), length(unique(seeds_present)))
})