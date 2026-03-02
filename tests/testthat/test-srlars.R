test_that("srlars runs correctly on synthetic data", {

    # --- Setup: Generate Simple Data ---
    set.seed(0)
    n <- 50
    p <- 20
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", 1:p) # Ensure names exist
    beta <- c(rep(2, 3), rep(0, p - 3))
    y <- as.numeric(x %*% beta + rnorm(n))

    # --- Test 1: Robust Execution (Default) ---
    # Should run without error
    fit <- srlars(x, y, n_models = 3, tolerance = 0.01, robust = TRUE)

    # Check Object Structure
    expect_s3_class(fit, "srlars")
    expect_equal(fit$n_models, 3)
    expect_true(fit$robust)
    expect_type(fit$active.sets, "list")
    expect_type(fit$coefficients, "list")

    # Check DDC object presence (needed for predict)
    expect_false(is.null(fit$ddc.object))

    # --- Test 2: Non-Robust Execution ---
    fit_nr <- srlars(x, y, n_models = 2, robust = FALSE)

    expect_s3_class(fit_nr, "srlars")
    expect_false(fit_nr$robust)
    expect_true(is.null(fit_nr$ddc.object))

    # --- Test 3: Coefficients Method ---
    # Should return vector of length p + 1 (intercept)
    coefs <- coef(fit)

    # FIX: Use expect_type for base types ("double" for numeric)
    expect_type(coefs, "double")
    expect_length(coefs, p + 1)

    # --- Test 4: Prediction Method ---
    # Test with new data
    newx <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
    colnames(newx) <- colnames(x) # Names must match for robust prediction

    # Robust prediction (Dynamic DDC)
    preds <- predict(fit, newx, dynamic = TRUE)
    expect_type(preds, "double")
    expect_length(preds, 10)

    # Non-robust prediction
    preds_nr <- predict(fit_nr, newx)
    expect_length(preds_nr, 10)

    # --- Test 5: Input Checks ---
    # Should error on mismatched dimensions
    expect_error(srlars(x, y[1:(n-1)]))
})
