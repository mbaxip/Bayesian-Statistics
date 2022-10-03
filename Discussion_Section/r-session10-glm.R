# MA 578 - Bayesian Statistics
# Bayesian GLMs, logistic regression

library(rstan)
library(bayesplot)

bioassay <- read.csv("data/bioassay.csv")
fit <- glm(cbind(deaths, n - deaths) ~ logdose, data = bioassay,
           family = binomial)
summary(fit)
# visualizing:
inv_logit <- plogis
x <- with(bioassay, seq(min(logdose), max(logdose), length = 100))
with(bioassay, plot(logdose, deaths / n, pch = 19))
lines(x, inv_logit(fit$coef[1] + fit$coef[2] * x), lwd = 2)


sm <- stan_model(model_code = "data {
  int<lower=0> N; int<lower=0> p;
  int<lower=0> n[N]; int<lower=0> y[N];
  matrix[N, p] X;
}
parameters { vector[p] beta; }
model { y ~ binomial_logit(n, X * beta); }")

y <- bioassay$deaths; n <- bioassay$n
X <- model.matrix(~ logdose, data = bioassay)
N <- nrow(X); p <- ncol(X)
sf <- sampling(sm, data = c("y", "n", "X", "N", "p"))
monitor(sf)
mcmc_trace(sf)
mcmc_acf(sf)
mcmc_dens_overlay(sf)

library(rstanarm)
sa <- stan_glm(cbind(deaths, n - deaths) ~ logdose, data = bioassay,
               family = binomial,
               prior_intercept = NULL, prior = NULL)
monitor(sa$stanfit) # equivalent

# visualizing
sims <- as.array(sf)
with(bioassay, plot(logdose, deaths / n, type = "n"))
for (it in 1:(dim(sims)[1])) {
  lines(x, inv_logit(sims[it, 1, 1] + sims[it, 1, 2] * x), col = "gray")
}
with(bioassay, points(logdose, deaths / n, pch = 19))
lines(x, inv_logit(fit$coef[1] + fit$coef[2] * x), lwd = 2)
