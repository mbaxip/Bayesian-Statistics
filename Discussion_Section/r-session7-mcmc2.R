# MA 578: R session, 11/7/19
# Bayesian computations: Markov chain Monte Carlo methods

library(rstan)
library(bayesplot)

# [ Normal hierarchical model (BDA Section 5.5) ]
# yij | theta_j, sigma^2 ~ind N(theta_j, sigma ^ 2)
# theta_j | mu, tau^2 ~iid N(mu, tau^2)
# P(mu, sigma^2, tau^2) propto 1 / sigma^2 * 1 / tau

# data and sufficient stats
diet <- read.csv("data/diet.csv", header = TRUE, comment = "#")
n <- with(diet, tapply(diet, diet, length))
J <- length(n)
idx <- c(0, cumsum(n)[-J])
N <- sum(n)
y <- diet$time

ns <- 2000 # number of mcmc samples
sm <- stan_model("diet.stan")
sf <- sampling(sm, iter = ns, chains = 4,
            data = c("J", "n", "idx", "N", "y"))

monitor(sf)
sims <- as.array(sf)
mcmc_trace(sims)
mcmc_acf(sims)
mcmc_dens_overlay(sims)
mcmc_hist(sims)


# [ Binomial hierarchical model (BDA Section 3.3) ]
# y_j | theta_j ~ind Binom(n_j, theta_j)
# theta_j | alpha, beta ~iid Beta(alpha, beta)
# P(alpha, beta) propto (alpha + beta)^(-2.5)

tumor <- read.csv("data/rattumor.csv", header = TRUE)
J <- nrow(tumor); n <- tumor$n; y <- tumor$y

ns <- 2000 # number of mcmc samples
sm <- stan_model("tumor.stan")
sf <- sampling(sm, iter = ns, chains = 4, data = c("J", "n", "y"))

monitor(sf)
sims <- as.array(sf)
params <- c(paste0("theta[", 1:3, "]"), "rho", "kappa", "lp__")
mcmc_trace(sims, params)
mcmc_acf(sims, params)
mcmc_dens_overlay(sims, params)
mcmc_hist(sims, params)
