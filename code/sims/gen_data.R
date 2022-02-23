## generate data for missing data feature selection simulation

expit <- function(x) exp(x) / (1 + exp(x))
## generate X data
#' @param n the sample size
#' @param p the number of covariates
#' @param rho the correlation
#' @param x_dist the distribution of x
#' @return matrix of covariates
gen_x <- function(n, p, rho = 0, x_dist = "normal") {
  if (x_dist == "normal") {
    sig <- matrix(rho, p, p)
    diag(sig) <- 1  
    x <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sig)
  } else if (x_dist == "nonnormal") {
    x1 <- rnorm(n, mean = 0.5, sd = 1)
    x2 <- rbinom(n, size = 1, prob = 0.5)
    x3 <- rweibull(n, shape = 1.75, scale = 1.9)
    log_x4 <- rnorm(n, mean = 0.5, sd = 0.5)
    x4 <- exp(log_x4)
    x5 <- rbinom(n, size = 1, prob = 0.5)
    x6 <- rnorm(n, mean = 0.25, sd = 1)
    if (p > 6) {
      sig <- matrix(rho, p - 6, p - 6)
      diag(sig) <- 1
      other_x <- MASS::mvrnorm(n = n, mu = rep(0, p - 6), Sigma = sig)
    } else {
      other_x <- NULL
    }
    x <- cbind(x1, x2, x3, x4, x5, x6, other_x)
  } else {
    stop("The entered distribution for x is not currently supported. Please enter one of 'normal' or 'nonnormal'.")
  }
  colnames(x) <- rep("", ncol(x))
  return(x)
}
## generate Y data
## @param x the covariates
## @param func the function that makes the linear predictor (e.g., function(x, beta) x%*%beta for a linear model)
## @param family specifies the link function ('gaussian' is identity link, 'binomial' is logit link)
gen_y <- function(x, func, family) {
  n <- dim(x)[1]
  linear_predictor <- func(x)
  if (family == "gaussian") {
    return(linear_predictor + rnorm(n, 0, 1))
  } else if (family == "binomial-logit" | family == "binomial") {
    return(rbinom(n, 1, prob = expit(linear_predictor)))
  } else if (family == "binomial-probit") {
    return(as.numeric((linear_predictor + rnorm(n, 0, 1)) > 0))  
  } else {
    stop("This function currently only works for family = 'gaussian' or family = 'binomial'. Please specify one of these families.")
  }
}
## generate a full dataset
gen_data <- function(n, p, rho, func, family, x_dist) {
  x <- gen_x(n = n, p = p, rho = rho, x_dist = x_dist)
  y <- gen_y(x = x, func = func, family = family)
  tib <- tibble::tibble(data.frame(y = y), as.data.frame(x))
  tib
}
## ----------------------------------------------------------------
## generate data from a linear outcome model
## ----------------------------------------------------------------
gen_data_lm <- function(n, p, rho, beta, family, x_dist) {
  return(gen_data(n = n, p = p, rho = rho, func = function(x) cbind(1, x)%*%beta, family = family, x_dist = x_dist))
}
## ----------------------------------------------------------------
## generate data from a nonlinear outcome model
## specifically, 10 features matter, and the rest don't
## sparse additive model
## ----------------------------------------------------------------
# center and scale based on population mean and variance
center_scale <- function(x, mean, sd) {
  (x - mean) / sd
}
f1 <- function(x) sin(pi / 4 * x)
f2 <- function(x, y) x * y
f3 <- function(x) tanh(x)
f4 <- function(x) cos(pi / 4 * x)
f5 <- function(x) (x ^ 2 + 1) ^ (-1) 
f6 <- function(x) (-1) * tanh(x)
nl_conditional_mean <- function(x, beta, x_dist) {
  centered_x <- apply(x, 2, function(col) center_scale(col, mean = 0, sd = 1))
  if (x_dist == "nonnormal") {
    centered_x[, 1] <- apply(x[, 1, drop = FALSE], 2, center_scale, mean = 0.5, sd = 1)
    centered_x[, 2] <- apply(x[, 2, drop = FALSE], 2, center_scale, mean = 0.5, sd = sqrt(0.5))
    centered_x[, 3] <- apply(x[, 3, drop = FALSE], 2, center_scale, 
                             mean = 1.9 * gamma(1 + 1 / 1.75), 
                             sd = sqrt(1.9 ^ 2 * (gamma(1 + 2 / 1.75) - gamma(1 + 1 / 1.75) ^ 2)))
    centered_x[, 4] <- apply(x[, 2, drop = FALSE], 2, center_scale, 
                             mean = exp(.5 + .5 ^ 2 / 2),
                             sd = sqrt((exp(.5 ^ 2) - 1) * exp(1 + .5 ^ 2)))
    centered_x[, 5] <- apply(x[, 5, drop = FALSE], 2, center_scale, mean = 0.5, sd = sqrt(0.5))
    centered_x[, 6] <- apply(x[, 6, drop = FALSE], 2, center_scale, mean = .25, sd = 1)
  } 
  new_x <- centered_x
  new_x[, 1] <- f1(centered_x[, 1])
  new_x[, 2] <- f2(centered_x[, 2], centered_x[, 3])
  new_x[, 3] <- f3(centered_x[, 3])
  new_x[, 4] <- f4(centered_x[, 4])
  new_x[, 5] <- f2(centered_x[, 5], centered_x[, 1])
  new_x[, 6] <- f6(centered_x[, 6])
  cbind(1, as.matrix(new_x)) %*% beta
}
gen_data_nlm <- function(n, p, rho, beta, family, x_dist) {
  return(gen_data(n = n, p = p, rho = rho, func = function(x) nl_conditional_mean(x, beta, x_dist), family = family, x_dist = x_dist))
}