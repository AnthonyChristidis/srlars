# srlars 1.0.0
* Initial stable release of package.

# srlars 1.0.1
* Fix error related to exhaustion of number of predictors.

# srlars 2.0.0
* Implementation of new algorithm for FSCRE.

# srlars 2.0.1
* Fix README to reflect new version of `srlars` function.

# srlars 3.0.0
* **Breaking change:** redesigned `srlars()` interface. Replaced `robust=TRUE/FALSE` with explicit preprocessing and CV controls (`x_preprocess`, `y_preprocess`, `cor_estimator`, `cv_preprocess`, `cv_fit`, `cv_loss`, `cv_folds`).
* **No target leakage by default:** predictor preprocessing is now performed on `X` only (no joint preprocessing of `[y, X]`).
* Added **wrapping-based robust, PSD correlations** (`cor_estimator = "wrap"`), based on `cellWise::wrap()`.
* Added **leakage-free foldwise CV preprocessing** option (`cv_preprocess = "foldwise"`) using `cellWise::DDCpredict()` and foldwise response transforms.
* Added robust arbiter options: **Huber/trimmed/MSE scoring** and optional **Huber IRLS fitting** inside the CV loop (`cv_fit = "huber"`).
* Updated internal selection loop stopping logic to require **strictly positive CV improvement** before accepting a variable.

# srlars 3.0.1
* Make `cv_fit = "huber"` the default.

# srlars 3.1.0
* Fixed `predict.srlars()`'s dynamic DDC-cleaning of new data: `object$robust` is now set by
  `srlars()` (it was previously always missing, so the DDC-cleaning branch was dead code and
  `dynamic = TRUE` had no effect), and the stale dummy-response-column augmentation before
  calling `cellWise::DDCpredict()` was removed to match the fact that `x_preprocess = "ddc"` is
  fit on the predictors alone.
* Added `max_share` argument to `srlars()`: the maximum number of sub-models (1 to `n_models`)
  a given variable may appear in. Default is `1`, reproducing the original fully-disjoint
  behavior exactly. When `1 < max_share < n_models`, each sub-model's first selected variable is
  forced to be distinct across sub-models, preventing several sub-models from redundantly
  duplicating the same strongest cold-start predictor; sharing is only permitted for variables
  added after a sub-model's first pick. This restriction is lifted entirely at
  `max_share = n_models` (sub-models are then free to become identical).
* Added `n_min` argument to `srlars()`: the minimum number of variables each sub-model is
  guaranteed to receive (subject to availability), even if a candidate doesn't clear the usual
  positive-benefit/`tolerance` requirement. Default is `NULL` (no floor, original behavior). Never
  bypasses the `max_share`/diversity pool restrictions -- only the CV-benefit acceptance
  requirement is relaxed for sub-models below the floor.
* Added `cv.srlars()`: chooses `max_share` by (outer) cross-validation and returns the ensemble
  refit at the cross-validated optimum. `coef()` and `predict()` work directly on the result via
  the existing `coef.srlars()`/`predict.srlars()` methods.