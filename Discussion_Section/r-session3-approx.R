# MA 578: R session, 10/01/19

# [ 1 ]
# y_i | theta ~ Po(x_i * theta), theta ~ Gamma(alpha, beta)
# Train example, as in lecture
x <- c(1, 2, 4)
y <- c(2, 4, 16)
alpha <- 1; beta <- 0 # "flat" prior
t_mode <- (alpha + sum(y) - 1) / (beta + sum(x))
t_var <- t_mode ^ 2 / (alpha + sum(y) - 1)

t <- seq(1, 6, length = 100)
plot(t, dgamma(t, alpha + sum(y), beta + sum(x)), type = "l",
     xlab = expression(theta), ylab = "density")
lines(t, dnorm(t, t_mode, sqrt(t_var)), lty = 2) # Laplace approx
abline(v = t_mode, lty = 3)


# [ 2 ]
# y_i | alpha, beta ~ Binomial(n_i, inv.logit(alpha + beta * x_i))
# P(alpha, beta) propto 1
# Extended analysis with
# (i) flat prior
# (ii) normal approximation

# read data (Table 3.1 in BDA)
bioassay <- read.csv("data/bioassay.csv", header = TRUE)

# check frequentist estimates
gb <- glm(cbind(deaths, n - deaths) ~ logdose,
          family = binomial, data = bioassay)
summary(gb) # check posterior mode and standard deviation to define grid range

m <- 100 # number of grid subdivisions
alpha <- seq(-2, 5, length = m) # (-5, 10) in Gelman
beta <- seq(-5, 25, length = m) # (-10, 40) in Gelman

# [ (i) flat prior ]
ilogit <- function (x) 1 / (1 + exp(-x))
lhood1 <- function (a, b) { # log-likelihood
  p <- ilogit(a + b * bioassay$logdose)
  sum(bioassay$deaths * log(p) + (bioassay$n - bioassay$deaths) * log(1 - p))
}

lj <- matrix(nrow = m, ncol = m) # log-joint grid
for (ia in 1:m)
  for (ib in 1:m)
    lj[ia, ib] <- lhood1(alpha[ia], beta[ib])
pab <- exp(lj - log(sum(exp(lj)))) # posterior on alpha, beta
# plot
image(alpha, beta, pab)
points(coef(gb)[1], coef(gb)[2], pch = 3)
contour(alpha, beta, pab, add = TRUE)


# [ (ii) normal approximation ]
precn <- solve(summary(gb)$cov.unscaled) # observed information
# here's another way to (re-)compute the information:
X <- cbind(rep(1, nrow(bioassay)), bioassay$logdose) # design matrix
# estimated proportions [check `predict(gb, type='response')`]:
p <- ilogit(drop(X %*% coef(gb)))
W <- bioassay$n * p * (1 - p) # variance weights
precn <- crossprod(X, W * X) # X' %*% diag(W) %*% X

lhoodn <- function (a, b) { # approximate log-likelihood
  v <- matrix(c(a, b) - coef(gb), ncol=1) # mean is mle, coef(gb)
  -.5 * crossprod(v, precn) %*% v
}

ln <- matrix(nrow = m, ncol = m) # log-lhood grid
for (ia in 1:m)
  for (ib in 1:m)
    ln[ia, ib] <- lhoodn(alpha[ia], beta[ib])
pabn <- exp(ln - log(sum(exp(ln)))) # approx. posterior on alpha, beta
# plot
image(alpha, beta, pab)
points(coef(gb)[1], coef(gb)[2], pch = 3)
contour(alpha, beta, pabn, add = TRUE) # compare approx contour with exact image

# marginal comparison
pa <- rowSums(pab) # exact marginal on alpha
pan <- rowSums(pabn) # approx marginal
ib_max <- which(pab == max(pab), arr.ind = TRUE)[, 2]
pam <- pab[, ib_max] / sum(pab[, ib_max]) # profile at beta^
plot(alpha, pam, type = "l", xlab = expression(alpha), ylab = "density")
lines(alpha, pa, col = "gray") # marginal on alpha
lines(alpha, pan, lty = 2) # note the lack of skewness
abline(v = coef(gb)[1], lty = 3) # posterior mode for alpha
