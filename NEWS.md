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