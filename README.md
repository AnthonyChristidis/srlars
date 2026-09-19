[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/srlars)](https://cran.r-project.org/package=srlars)
[![CRAN Data](https://www.r-pkg.org/badges/last-release/srlars)](https://cran.r-project.org/package=srlars)
[![Downloads](https://cranlogs.r-pkg.org/badges/srlars)](https://cran.r-project.org/package=srlars)
[![arXiv](https://img.shields.io/badge/arXiv-2603.20940-b31b1b.svg)](https://arxiv.org/abs/2603.20940)

# srlars: Fast and Scalable Cellwise-Robust Ensembles

High-dimensional data are often affected by **cellwise contamination**: individual
*cells* of the predictor matrix deviate from the underlying structure without
necessarily making the whole observation an outlier. Even a small fraction of
contaminated cells can propagate across many observations, which is enough to
mislead both classical variable selection methods and robust methods designed only
for casewise (whole-observation) outliers.

`srlars` implements the **Fast and Scalable Cellwise-Robust Ensemble (FSCRE)**
algorithm: a competitive ensemble of sparse sub-models built on a cellwise-robust
foundation (Detect Deviating Cells imputation and wrapping-based robust
correlations), constructed via a robust Least-Angle-Regression proposer and a
cross-validation arbiter, then refit with robust MM-estimators. The method and its
theoretical properties are described in:

> Christidis, A., Pyneeandee, J., and Cohen Freue, G. (2026). *Fast and Scalable
> Cellwise-Robust Ensembles for High-Dimensional Data*.
> [arXiv:2603.20940](https://arxiv.org/abs/2603.20940)

## Key features

- **Cellwise-robust foundation** -- predictors are cleaned with `cellWise::DDC()`
  and correlations are estimated with the wrapping transform, so estimation stays
  reliable when individual cells (not whole rows) are contaminated.
- **Competitive ensemble construction** -- `n_models` sub-models compete for
  variables each round via cross-validated predictive improvement, rather than
  being built independently.
- **Controllable variable sharing** -- `max_share` sets how many sub-models a given
  variable may appear in, from fully disjoint sub-models (the default) to
  unrestricted sharing.
- **Minimum sub-model size** -- `n_min` guarantees each sub-model reaches a minimum
  number of variables even when the ensemble-wide stopping rule would otherwise cut
  it short.
- **Automatic tuning** -- `cv.srlars()` chooses `max_share` by cross-validation on
  held-out ensemble prediction error, instead of comparing values by hand.

## Installation

You can install the **stable** version from [CRAN](https://cran.r-project.org/package=srlars):

```r
install.packages("srlars", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/srlars):

```r
library(devtools)
devtools::install_github("AnthonyChristidis/srlars")
```

## Quick start

```r
library(srlars)

# x, y: a (possibly cellwise-contaminated) high-dimensional training set
fit <- srlars(x, y,
             n_models = 5,       # ensemble size
             x_preprocess = "ddc",
             y_preprocess = "wrap",
             cor_estimator = "wrap",
             cv_fit = "huber",
             cv_loss = "huber")

coef(fit)          # ensemble-averaged coefficients
predict(fit, newx) # ensemble-averaged predictions

# Choose max_share automatically instead of setting it by hand:
cv_fit <- cv.srlars(x, y, n_models = 5)
coef(cv_fit) # coef()/predict() work directly on the cross-validated fit
```

For a complete walkthrough -- simulating cellwise-contaminated data, fitting
`srlars()`, and comparing `max_share`/`n_min` on the same dataset -- see the
package vignette:

```r
vignette("srlars", package = "srlars")
```

## License

This package is free and open source software, licensed under GPL (&gt;= 2).
