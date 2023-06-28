[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/srlars)](https://cran.r-project.org/package=srlars)
[![CRAN Data](https://www.r-pkg.org/badges/last-release/srlars)](https://cran.r-project.org/package=srlars) 
[![Downloads](http://cranlogs.r-pkg.org/badges/srlars)](https://cran.r-project.org/package=srlars)

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
# Simulation parameters
n <- 50
p <- 500
rho <- 0.5
rho.inactive <- 0.2
group.size <- 25
p.active <- 100
snr <- 1
contamination.prop <- 0.2

# Setting the seed
set.seed(0)

# Block Correlation
sigma.mat <- matrix(0, p, p)
sigma.mat[1:p.active, 1:p.active] <- rho.inactive
for(group in 0:(p.active/group.size - 1))
  sigma.mat[(group*group.size+1):(group*group.size+group.size),(group*group.size+1):(group*group.size+group.size)] <- rho
diag(sigma.mat) <- 1

# Simulation of beta vector
true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))

# Setting the SD of the variance
sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))

# Simulation of test data
m <- 2e3
x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)

# Simulation of uncontaminated data
x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
y <- x %*% true.beta + rnorm(n, 0, sigma)

# Contamination of data
contamination_indices <- 1:floor(n*contamination.prop)
k_lev <- 2
k_slo <- 100
x_train <- x
y_train <- y
beta_cont <- true.beta
beta_cont[true.beta!=0] <- beta_cont[true.beta!=0]*(1 + k_slo)
beta_cont[true.beta==0] <- k_slo*max(abs(true.beta))
for(cont_id in contamination_indices){

  a <- runif(p, min = -1, max = 1)
  a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
  x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
  y_train[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
}

# srlars models
srlars_fit <- srlars(x_train, y_train,
                     n_models = 5,
                     model_saturation = c("fixed", "p-value")[1],
                     alpha = 0.05, model_size = n-1,
                     robust = TRUE,
                     compute_coef = TRUE,
                     en_alpha = 1/4)
srlars_preds <- predict(srlars_fit, newx = x_test,
                        group_index = 1:srlars_fit$n_models,
                        dynamic = FALSE)
srlars_coefs <- coef(srlars_fit, group_index = 1:srlars_fit$n_models)
sens_srlars <- sum(which((srlars_coefs[-1]!=0)) <= p.active)/p.active
spec_srlars <- sum(which((srlars_coefs[-1]!=0)) <= p.active)/sum(srlars_coefs[-1]!=0)
mspe_srlars <- mean((y_test - srlars_preds)^2)/sigma^2
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
