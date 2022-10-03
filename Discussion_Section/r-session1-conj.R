# MA 578: R session, 09/19/19
# y_i | theta ~ Poisson(theta * x_i)

# [ Prior Elicitation ]
# assume conjugacy: theta ~ Gamma(alpha, beta)
# "on average, 4 cars per hour" -> E[theta] = alpha / beta = 4
# "most likely between 1 to 6 cars per hour" -> P(1 <= theta <= 6) = .95
# find beta numerically by finding an approximate root to 'f'
f <- function (beta)
  pgamma(6, 4 * beta, beta) - pgamma(1, 4 * beta, beta) - .95
beta <- round(uniroot(f, lower = 0, upper = 10)$root, 0)
b <- seq(0, 10, length = 100)
plot(b, f(b), type = "l")
abline(h = 0, lty = 2); abline(v = beta)
alpha <- 4 * beta

# [ Bayesian Analysis ]
# data:
x <- c(1, 2, 4); y <- c(2, 4, 16)
# check posterior graphically
t <- seq(0, 8, length = 100) # theta
plot(t, dgamma(t, alpha + sum(y), beta + sum(x)),
     type = "l", xlab = expression(theta), ylab = "density")
# 95% credible interval
qgamma(c(.025, .975), alpha + sum(y), beta + sum(x))

# describing posterior by sampling:
ns <- 1000 # #samples
ts <- rgamma(ns, alpha + sum(y), beta + sum(x))
hist(ts, prob = TRUE)
lines(t, dgamma(t, alpha + sum(y), beta + sum(x)))
mean(ts)
quantile(ts, c(.025, .975))

# compare to prior
plot(t, dgamma(t, alpha + sum(y), beta + sum(x)),
     type = "l", xlab = expression(theta), ylab = "density")
lines(t, dgamma(t, alpha, beta), lty = 2) # prior
lines(t, dgamma(t, sum(y) + 1, sum(x)), lty = 3) # likelihood

# posterior predictive
help(dnbinom) # check negative binomial distribution
xt <- 1 # x tilde
yt <- 0:10
plot(yt,
     dnbinom(yt, alpha + sum(y), 1 - xt / (xt + beta + sum(x))),
     type = "h", xlab = expression(tilde(y)),
     ylab = expression(P(tilde(y)~"|"~y)))
# 95% posterior predictive interval
qnbinom(c(.025, .975), alpha + sum(y), 1 - xt / (xt + beta + sum(x)))

# using sampling:
ys <- rnbinom(ns, alpha + sum(y), 1 - xt / (xt + beta + sum(x)))
mean(ys)
quantile(ys, c(.025, .975))
# conditional sampling: sample theta | y, then y.tilde | theta
ys <- rpois(ns, rgamma(ns, alpha + sum(y), beta + sum(x)))
# or, in two steps,
ts <- rgamma(ns, alpha + sum(y), beta + sum(x))
ys <- rpois(ns, xt * ts)
mean(ys) # compare to marginal sample from negative binomial
quantile(ys, c(.025, .975))
