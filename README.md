[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/srlars)](https://cran.r-project.org/package=srlars)
[![CRAN Data](https://www.r-pkg.org/badges/last-release/srlars)](https://cran.r-project.org/package=srlars) 
[![Downloads](https://cranlogs.r-pkg.org/badges/srlars)](https://cran.r-project.org/package=srlars)

srlars
================

This package provides functions for performing split robust least angle regression.

---------------------------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=srlars).

```{r installation, eval = FALSE}
install.packages("srlars", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/srlars)

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/srlars")
```

### Usage

``` r
library(srlars)
library(mvnfast)

# --- 1. Simulation Parameters ---

n <- 50
p <- 100
rho.within <- 0.8
rho.between <- 0.2
p.active <- 20
group.size <- 5
snr <- 3
contamination.prop <- 0.1

# Setting the seed
set.seed(0)

# --- 2. Data Generation ---

# Block correlation structure
sigma.mat <- matrix(0, p, p)
sigma.mat[1:p.active, 1:p.active] <- rho.between
for(group in 0:(p.active/group.size - 1))
  sigma.mat[(group*group.size+1):(group*group.size+group.size),
  (group*group.size+1):(group*group.size+group.size)] <- rho.within
diag(sigma.mat) <- 1

# True coefficient vector
true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))

# Noise level
sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))

# Generate uncontaminated training data
x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
y <- x %*% true.beta + rnorm(n, 0, sigma)

# Generate test data
m <- 2e3
x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)

# --- 3. Introduce Contamination ---

# Cellwise contamination
contamination_indices <- sample(1:(n * p), round(n * p * contamination.prop))
x_train <- x
x_train[contamination_indices] <- runif(length(contamination_indices), -10, 10)
y_train <- y

# --- 4. Fit srlars Model ---

# Fit the FSCRE ensemble
# We use 5 sub-models and robust initialization
fit <- srlars(x_train, y_train,
              n_models = 5,
              tolerance = 0.01,
              robust = TRUE,
              compute_coef = TRUE)

# --- 5. Prediction and Evaluation ---

# Predict on new data
# By default, this uses dynamic robust imputation for the new data
preds <- predict(fit, newx = x_test)

# Evaluate MSPE
mspe <- mean((y_test - preds)^2) / sigma^2
print(paste("MSPE:", round(mspe, 3)))

# Extract Coefficients (averaged over the ensemble)
coefs <- coef(fit)

# Variable Selection Metrics
selected_indices <- which(coefs[-1] != 0)
true_indices <- which(true.beta != 0)

recall <- length(intersect(selected_indices, true_indices)) / length(true_indices)
precision <- length(intersect(selected_indices, true_indices)) / length(selected_indices)

print(paste("Recall:", round(recall, 3)))
print(paste("Precision:", round(precision, 3)))
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
