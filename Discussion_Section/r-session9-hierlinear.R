# MA 578 - Bayesian Statistics
# Bayesian Hierarchical Linear Models

library(rstan)
library(bayesplot)
source("~/Documents/BU_2019_Fall/bslm.R")

# [ Setup ]
stroke <- within(read.csv("~/Documents/BU_2019_Fall/stroke.csv"), subject <- factor(subject))

library(lattice)
xyplot(score ~ week | group, data = stroke)
xyplot(score ~ week | subject, data = stroke)
xyplot(score ~ week | subject, data = stroke,
       panel = function(...) {
         panel.xyplot(type = "p", ...)
         panel.lmline(lty = 2, ...)
       })


# setup design matrices
y <- stroke$score
X <- model.matrix(~ -1 + subject + subject:week, data = stroke)
J <- length(levels(stroke$subject))
n <- nrow(X)
p <- ncol(X) / J
G <- length(levels(stroke$group))
stroke_j <- aggregate(. ~ subject + group, data = stroke, length)
stroke_subject <- data.frame(group = stroke_j$group,
    coef = c(rep("_intercept", J), rep("_week", J)))
X_beta <- model.matrix(~ -1 + coef + coef:group, data = stroke_subject)

# priors
p_alpha <- ncol(X_beta); p_beta <- J * p
alpha_prior <- list(mean = rep(0, p_alpha), precision = rep(0, p_alpha)) # "fixed"
sigma2_prior <- list(df = 0, scale = 0)
tau2_prior <- list(df = 0, scale = 0)

# [ Hierarchical model ]
iter <- 2000; warmup <- floor(iter / 2)
nchains <- 4
params <- c(colnames(X_beta), "sigma", "tau", "lp__")
sims <- mcmc_array(iter - warmup, nchains, params)

bl <- lm(y ~ X - 1)
beta <- coef(bl)
sigma <- sqrt(deviance(bl) / n)
al <- lm(beta ~ X_beta - 1)
tau <- sqrt(deviance(al) / p_beta)
for (chain in 1:nchains) {
  for (it in 1:iter) {
    alpha <- bslm_sample1(beta, X_beta, tau, alpha_prior)
    beta <- bslm_sample1(y, X, sigma,
                         list(mean = X_beta %*% alpha, precision = 1 / tau ^ 2))
    sigma <- sqrt(rinvchisq(1, sigma2_prior$df + n,
                             sigma2_prior$df * sigma2_prior$scale +
                               crossprod(y - X %*% beta)))
    tau <- sqrt(rinvchisq(1, tau2_prior$df + p_beta,
                          tau2_prior$df * tau2_prior$scale +
                            crossprod(beta - X_beta %*% alpha)))
    target <- sum(dnorm(y, X %*% beta, sigma, log = TRUE)) +
      sum(dnorm(beta, X_beta %*% alpha, tau, log = TRUE)) -
      log(sigma) - log(tau)

    if (it > warmup)
      sims[it - warmup, chain, ] <- c(alpha, sigma, tau, target)
  }
}

monitor(sims)
mcmc_trace(sims)
mcmc_acf(sims)
mcmc_dens_overlay(sims)

# or, in Stan:
sm <- stan_model(model_code = "data {
  int<lower=0> n; int<lower=0> J; int<lower=0> p; int<lower=0> G;
  vector[n] y; matrix[n, p * J] X; matrix[p * J, p * G] X_beta;
}
parameters {
  vector[p * J] beta; vector[p * G] alpha;
  real<lower=0> sigma; real<lower=0> tau;
}
model {
  y ~ normal(X * beta, sigma);
  beta ~ normal(X_beta * alpha, tau);
  target += -log(sigma) - log(tau);
}")
sf <- sampling(sm, data = c("n", "J", "p", "G", "y", "X", "X_beta"))
monitor(sf)


# [ model assessment ]
W <- with(stroke, data.frame(split(score, week)))
corW <- cor(W)
nweek <- nrow(corW)
indW <- which(outer(1:nweek, 1:nweek, `-`) > 0, arr.ind = TRUE) # lower triangle
lag <- indW[,1] - indW[,2]
lcorW <- log(abs(cor(W))[indW])
rho <- coef(lm(lcorW ~ -1 + lag))
# simple autocorrelation not a good fit, but can work as a test function:
boxplot(lcorW ~ lag)
points(1:(nweek - 1), rho * (1:(nweek-1)), pch = 19, col = "red")

# log(estimated autocorrelation)
cor_test <- function (s) { # 's' are scores
  W <- data.frame(split(s, stroke$week))
  lcorW <- log(abs(cor(W))[indW])
  coef(lm(lcorW ~ -1 + lag))
}

nsamples <- iter - warmup
Z <- X %*% X_beta
sims0 <- bslm_sample(y, Z) # naive model on fixed effects only
ct0 <- numeric(nsamples)
for (it in 1:nsamples) {
  alpha <- sims0[it, 1, 1:p_alpha]
  sigma <- sims0[it, 1, p_alpha + 1]
  s <- rnorm(n, Z %*% alpha, sigma)
  ct0[it] <- cor_test(s)
}

ct <- numeric(nsamples)
for (it in 1:nsamples) {
  beta <- rnorm(p_beta, X_beta %*% sims[it, 1, 1:p_alpha], sims[it, 1, p_alpha + 2])
  s <- rnorm(n, X %*% beta, sims[it, 1, p_alpha + 1])
  ct[it] <- cor_test(s)
}

# hierarchical model captures auto-correlation better:
boxplot(ct0, ct, col = "gray"); abline(h = rho, col = "red")

# but it has very large variance on the slope
library(beanplot)
alpha <- sims[nsamples, 1, 1:p_alpha]
sigma <- sims[nsamples, 1, p_alpha + 1]
tau <- sims[nsamples, 1, p_alpha + 2]
beta <- rnorm(p_beta, X_beta %*% alpha, tau)
s <- rnorm(n, X %*% beta, sigma)
test <- data.frame(score = c(stroke$score, s),
                   group = rep(stroke$group, 2),
                   type = c(rep("observed", n), rep("replicated", n)))
beanplot(score ~ interaction(type, group), side = "both", data = test)
abline(h = 0, col = "red", lty = 2)

# comparison
library(lme4)
# E[y_{ijk}] = (alpha + beta_j + gamma * w_k + theta_j * w_k)
#               + (phi_i + psi_i * w_k)
re <- lmer(score ~ group * week + (1 + week | subject), data = stroke)
