# MA 578, Bayesian Statistics
# Bayesian Linear Models
# y | beta, sigma2 ~ N(X * beta, sigma2 * I_n)
# beta ~ N(beta0, S_inv0 ^ (-1))
# sigma2 ~ Inv-scaled-Chisq(nu, tau2)

# [ Auxiliary ]
# sample from inverse (scaled) chi-square with parameters `nu` and `tau2`;
# nu_tau2 = nu * tau2 for convenience
rinvchisq <- function (ns, nu, nu_tau2) 1 / rgamma(ns, nu / 2, nu_tau2 / 2)

mcmc_array <- function (ns, nchains = 1, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}


# [ Simple interface ]
bslm_fit <- function (y, x, prior_coef = NULL, prior_disp = NULL,
                      maxit = 25, epsilon = 1e-8) {
  nvars <- ncol(x); nobs <- nrow(x)
  dn <- colnames(x); if (is.null(dn)) dn <- paste0("x", 1L:nvars)
  if (is.null(prior_coef))
    prior_coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
  if (is.null(prior_disp))
    prior_disp <- list(df = 0, scale = 0)
  S_inv0 <- prior_coef$precision
  beta0 <- prior_coef$mean
  beta0 <- if (is.vector(S_inv0)) S_inv0 * beta0 else S_inv0 %*% beta0
  nu <- prior_disp$df
  nu_tau2 <- nu * prior_disp$scale

  rss <- sum((y - mean(y)) ^ 2) # intercept only
  sigma2 <- (nu_tau2 + rss) / (nu + nobs)
  for (iter in 1:maxit) {
    z <- crossprod(x, y) / sigma2 + beta0
    V_inv <- crossprod(x) / sigma2
    if (is.vector(S_inv0)) # diagonal precision?
      diag(V_inv) <- diag(V_inv) + S_inv0
    else
      V_inv <- V_inv + S_inv0
    C <- chol(V_inv)
    u <- backsolve(C, z, transpose = TRUE)
    coef <- drop(backsolve(C, u))
    rss_new <- sum((y - drop(x %*% coef)) ^ 2)
    sigma2 <- (nu_tau2 + rss_new) / (nu + nobs)

    rel_error <- abs((rss_new - rss) / rss)
    if (!is.infinite(rss_new) && (rel_error < epsilon)) # converged?
      break
    rss <- rss_new
  }
  names(coef) <- dn
  list(coef = coef, sigma = sqrt(sigma2), C = C)
}


bslm_sample <- function (y, x, prior_coef = NULL, prior_disp = NULL,
                         chains = 4, iter = 2000, warmup = floor(iter / 2)) {
  nvars <- ncol(x); nobs <- nrow(x)
  dn <- colnames(x); if (is.null(dn)) dn <- paste0("x", 1L:nvars)
  if (is.null(prior_coef))
    prior_coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
  if (is.null(prior_disp))
    prior_disp <- list(df = 0, scale = 0)
  S_inv0 <- prior_coef$precision
  beta0 <- prior_coef$mean
  beta0 <- if (is.vector(S_inv0)) S_inv0 * beta0 else S_inv0 %*% beta0
  nu <- prior_disp$df
  nu_tau2 <- nu * prior_disp$scale

  rss <- sum((y - mean(y)) ^ 2)
  sigma2 <- (nu_tau2 + rss) / (nu + nobs)
  sims <- mcmc_array(iter - warmup, chains, c(dn, "sigma", "lp__"))
  for (chain in 1:chains) {
    for (it in 1:iter) {
      z <- crossprod(x, y) / sigma2 + beta0
      V_inv <- crossprod(x) / sigma2
      if (is.vector(S_inv0)) # diagonal precision?
        diag(V_inv) <- diag(V_inv) + S_inv0
      else
        V_inv <- V_inv + S_inv0
      C <- chol(V_inv)
      u <- backsolve(C, z, transpose = TRUE)
      coef <- drop(backsolve(C, u + rnorm(nvars)))
      rss <- sum((y - drop(x %*% coef)) ^ 2)
      sigma2 <- rinvchisq(1, nu + nobs, nu_tau2 + rss)
      lp <- -((nu + nobs) / 2 + 1) * log(sigma2) - .5 * (nu_tau2 + rss) / sigma2
      if (it > warmup)
        sims[it - warmup, chain, ] <- c(coef, sqrt(sigma2), lp)
    }
  }
  sims
}


# y | beta, sigma2 ~ N(x * beta, sigma2 * I_n), beta ~ N(beta0, S_inv0^{-1}),
# sample from beta | sigma2, y ~ N(B * b, B) with
# B^{-1} = x' * x / sigma2 + S_inv0 and b = x' * y / sigma2 + S_inv0 * beta0
bslm_sample1 <- function (y, x, sigma, prior_coef) {
  beta0 <- prior_coef$mean; S_inv0 <- prior_coef$precision
  beta0 <- if (is.vector(S_inv0)) S_inv0 * beta0 else S_inv0 %*% beta0
  z <- crossprod(x, y) / sigma ^ 2 + beta0
  V_inv <- crossprod(x) / sigma ^ 2
  if (is.vector(S_inv0)) # diagonal precision?
    diag(V_inv) <- diag(V_inv) + S_inv0
  else
    V_inv <- V_inv + S_inv0
  C <- chol(V_inv)
  u <- backsolve(C, z, transpose = TRUE)
  drop(backsolve(C, u + rnorm(ncol(x))))
}
