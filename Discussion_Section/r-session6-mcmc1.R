# MA 578: R session, 10/31/19
# Bayesian computations: Markov chain Monte Carlo methods

# [ Auxiliary functions ]

# sample from inverse (scaled) chi-square with parameters `nu` and `tau2`;
# nu_tau2 = nu * tau2 for convenience
rinvsquare <- function (ns, nu, nu_tau2) 1 / rgamma(ns, nu / 2, nu_tau2 / 2)

mcmc_array <- function (ns, nchains = 1, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

plot_trace <- function (x, ...)
  plot(x, type = "l", col = "gray", xlab = "iteration", ...)
plot_hist <- function (x, ...)
  hist(x, col = "gray", border = "white", main = "", ...)


# [ Normal hierarchical model (BDA Section 5.5) ]
# yij | theta_j, sigma^2 ~ind N(theta_j, sigma ^ 2)
# theta_j | mu, tau^2 ~iid N(mu, tau^2)
# P(mu, sigma^2, tau^2) propto 1 / sigma^2 * 1 / tau

# data and sufficient stats
diet <- read.csv("data/diet.csv", header = TRUE, comment = "#")
nj <- with(diet, tapply(diet, diet, length))
yj <- with(diet, tapply(time, diet, mean))
s2j <- with(diet, tapply(time, diet, var))
S <- sum(s2j * (nj - 1))
J <- length(nj)
n <- sum(nj)

ns <- 10000 # number of mcmc samples
params <- c(paste0("theta", 1:J), "mu", "sigma2", "tau2", "lp__") 
sims <- mcmc_array(ns, params = params)

# init
theta <- with(diet, tapply(time, diet, mean))
mu <- mean(theta)
# iterate
for (is in 1:ns) {
  # [ sigma2 | theta, mu, tau2, y ]
  sigma2 <- rinvsquare(1, n, S + sum(nj * (yj - theta) ^ 2))
  # [ tau2 | theta, mu, sigma2, y ]
  tau2 <- rinvsquare(1, J - 1, sum((theta - mu) ^ 2))
  # [ theta | mu, sigma2, tau2, y ]
  theta <- rnorm(J, (nj * yj / sigma2 + mu / tau2) / (nj / sigma2 + 1 / tau2),
                 1 / sqrt(nj / sigma2 + 1 / tau2))
  # [ mu | theta, sigma2, tau2, y ]
  mu <- rnorm(1, mean(theta), sqrt(tau2 / J))
  # [ log posterior ]
  lp <- with(diet, sum(dnorm(time, theta[diet], sqrt(sigma2), log = TRUE))) +
    sum(dnorm(theta, mu, sqrt(tau2), log = TRUE)) - log(sigma2) - .5 * log(tau2)

  sims[is, 1, ] <- c(theta, mu, sigma2, tau2, lp)
}


r <- 8000:10000 # discard burn-in
op <- par(mfrow = c(2, 4))
plot_trace(sims[r, , J + 1], ylab = expression(mu))
plot_trace(sqrt(sims[r, , J + 2]), ylab = expression(sigma))
plot_trace(sqrt(sims[r, , J + 3]), ylab = expression(tau))
plot_trace(sims[r, , J + 4], ylab = "log-posterior")
for (j in 1:J)
  plot_trace(sims[r, , j], ylab = substitute(theta[j], list(j = j)))
par(op)

op <- par(mfrow = c(2, 4))
acf(sims[r, , J + 1], main = expression(mu))
acf(sqrt(sims[r, , J + 2]), main = expression(sigma))
acf(sqrt(sims[r, , J + 3]), main = expression(tau))
acf(sims[r, , J + 4], main = "log-posterior")
for (j in 1:J)
  acf(sims[r, , j], main = substitute(theta[j], list(j = j)))
par(op)

op <- par(mfrow = c(2, 4))
plot_hist(sims[r, , J + 1], xlab = expression(mu)); abline(v = mean(yj))
plot_hist(sqrt(sims[r, , J + 2]), xlab = expression(sigma))
plot_hist(sqrt(sims[r, , J + 3]), xlab = expression(tau))
plot_hist(sims[r, , J + 4], xlab = "log-posterior")
for (j in 1:J) {
  plot_hist(sims[r, , j], xlab = substitute(theta[j], list(j = j)))
  abline(v = yj[j])
}
par(op)


# [ Binomial hierarchical model (BDA Section 3.3) ]
# y_j | theta_j ~ind Binom(n_j, theta_j)
# theta_j | alpha, beta ~iid Beta(alpha, beta)
# P(alpha, beta) propto (alpha + beta)^(-2.5)

tumor <- read.csv("data/rattumor.csv", header = TRUE)
J <- nrow(tumor)
sigma_a <- .5 # randow walk sigma for alpha
sigma_b <- 4 # random walk sigma for beta

lprior <- function (a, b) -2.5 * log(a + b)
lhood_ab <- function (a, b)
  with(tumor,
    sum(lgamma(a + y) - lgamma(a) + lgamma(b + n - y) - lgamma(b) -
        (lgamma(a + b + n) - lgamma(a + b))))

ns <- 10000 # number of mcmc samples
params <- c(paste0("theta", 1:J), "alpha", "beta", "lp__") 
sims <- mcmc_array(ns, params = params)

alpha <- beta <- 1
for (is in 1:ns) {
  # [ theta | alpha, beta, y ]
  theta <- with(tumor, rbeta(J, alpha + y, beta + n - y))

  # [ alpha | theta, beta, y: randow walk ]
  alpha_c <- rnorm(1, alpha, sigma_a)
  if (alpha_c > 0) {
    log_r <- (alpha_c - alpha) * sum(log(theta)) +
              J * (lgamma(alpha_c + beta) - lgamma(alpha_c) -
                   (lgamma(alpha + beta) - lgamma(alpha))) +
              lprior(alpha_c, beta) - lprior(alpha, beta)
    if (log_r >= 0 || log(runif(1)) <= log_r) # accept?
      alpha <- alpha_c
  }

  # [ beta | theta, alpha, y: randow walk ]
  beta_c <- rnorm(1, beta, sigma_b)
  if (beta_c > 0) {
    log_r <- (beta_c - beta) * sum(log(1 - theta)) +
              J * (lgamma(alpha + beta_c) - lgamma(beta_c) -
                   (lgamma(alpha + beta) - lgamma(beta))) +
              lprior(alpha, beta_c) - lprior(alpha, beta)
    if (log_r >= 0 || log(runif(1)) <= log_r) # accept?
      beta <- beta_c
  }

  # [ log posterior ]
  lp <- with(tumor, sum(dbinom(y, n, theta, log = TRUE))) +
    sum(dbeta(theta, alpha, beta, log = TRUE)) + lprior(alpha, beta)

  sims[is, 1, ] <- c(theta, alpha, beta, lp)
}

r <- 8000:10000
op <- par(mfrow = c(3, 3))
plot_trace(sims[r, , J + 1], ylab = expression(alpha))
plot_trace(sims[r, , J + 2], ylab = expression(beta))
plot_trace(sims[r, , J + 3], ylab = "log-posterior")
acf(sims[r, , J + 1], main = expression(alpha))
acf(sims[r, , J + 2], main = expression(beta))
acf(sims[r, , J + 3], main = "log-posterior")
plot_hist(sims[r, , J + 1], xlab = expression(alpha))
plot_hist(sims[r, , J + 2], xlab = expression(beta))
plot_hist(sims[r, , J + 3], xlab = "log-posterior")
par(op)


# [ verify: numerical approximation to P(alpha, beta | y) ]
rho <- with(tumor, sum(y) / sum(n))
m <- 100
k <- seq(10, 20, length = m)
ks <- k[which.max(sapply(k, function (x) lhood_ab(rho * x, (1 - rho) * x)))]
am <- rho * ks; bm <- (1 - rho) * ks
a <- seq(0, 2 * am, length = m + 1)[-1]
b <- seq(0, 2 * bm, length = m + 1)[-1]

lab <- matrix(nrow = m, ncol = m)
for (ia in 1:m)
  for (ib in 1:m)
    lab[ia, ib] <- lprior(a[ia], b[ib]) + lhood_ab(a[ia], b[ib])
pab <- exp(lab - log(sum(exp(lab))))

image(a, b, pab); lines(sims[r, , J + 1], sims[r, , J + 2])
pa <- rowSums(pab); pb <- colSums(pab)
qqplot(sample(a, length(r), replace = TRUE, prob = pa), sims[r, , J + 1])
qqplot(sample(b, length(r), replace = TRUE, prob = pb), sims[r, , J + 2])

op <- par(mfrow = c(3, 3))
ind_j <- sample.int(J, 3)
for (j in ind_j) {
  label <- substitute(theta[j], list(j = j))
  plot_trace(sims[r, , j], ylab = label)
  acf(sims[r, , j], main = label)
  plot_hist(sims[r, , j], xlab = label)
  abline(v = with(tumor, y[j] / n[j]))
}
par(op)

