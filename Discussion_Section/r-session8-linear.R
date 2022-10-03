# MA 578 - Bayesian Statistics
# Bayesian Linear Models

library(rstan)
library(bayesplot)
library(beanplot)
source("bslm.R")
plot_hist <- function (x, ...)
  hist(x, col = "gray", border = "white", ...)


data(cars)
y <- cars$dist
X <- model.matrix(~ speed, data = cars)
n <- length(y); p <- ncol(X)

fit <- bslm_fit(y, X)
data.frame(coef = fit$coef, se = sqrt(diag(chol2inv(fit$C))))
summary(lm(dist ~ speed, data = cars))

sims <- bslm_sample(y, X)
monitor(sims)
mcmc_trace(sims)
mcmc_acf(sims)
mcmc_dens_overlay(sims)
mcmc_hist(sims)

# visualize
with(cars, plot(speed, dist, type = "n"))
ns <- dim(sims)[1]
for (s in 1:ns) abline(sims[s, 1, 1], sims[s, 1, 2], col = "gray")
with(cars, points(speed, dist, pch = 19))
abline(fit$coef[1], fit$coef[2], lwd = 2)

# compare to Stan
sm <- stan_model(model_code = "data {int<lower=0> n; int<lower=0> p; matrix[n, p] X; vector[n] y; }
                 parameters { vector[p] beta; real<lower=0> sigma; }
                 model { y ~ normal(X * beta, sigma); target += -log(sigma); }")
sf <- sampling(sm, data = c("n", "p", "X", "y"))


# model checking 
y_rep <- matrix(nrow = ns, ncol = n)
for (is in 1:ns)
  y_rep[is,] <- rnorm(n, X %*% sims[is, 1, 1:p], sims[is, 1, p + 1])

# visual check
boxplot(y_rep, outline = F); points(y, pch = 19, col = "red")
# test function: max(cook's distance)
hat <- diag(crossprod(backsolve(fit$C, t(X), transpose = T))) / fit$sigma ^ 2
ch <- hat / (1 - hat) ^ 2 / p
cook_dist <- function (std_res) ch * std_res ^ 2
cd0 <- numeric(ns); cd <- numeric(ns)
for (is in 1:ns) {
  r0 <- (y - X %*% sims[is, 1, 1:p]) / sims[is, 1, p + 1]
  cd0[is] <- max(cook_dist(r0))
  r <- (y_rep[is,] - X %*% sims[is, 1, 1:p]) / sims[is, 1, p + 1]
  cd[is] <- max(cook_dist(r))
}
mean(cd > cd0) # Bayesian "p-value"
beanplot(cd0, cd, side = "both")


# informative prior
pc <- list(mean = c(0, 0), precision = c(1, 0)) # shrink intercept to zero
fit <- bslm_fit(y, X, prior_coef = pc)
data.frame(coef = fit$coef, se = sqrt(diag(chol2inv(fit$C))))
sims <- bslm_sample(y, X, prior_coef = pc)
monitor(sims)
# visualize
with(cars, plot(speed, dist, type = "n"))
ns <- dim(sims)[1]
for (s in 1:ns) abline(sims[s, 1, 1], sims[s, 1, 2], col = "gray")
with(cars, points(speed, dist, pch = 19))
abline(fit$coef[1], fit$coef[2], lwd = 2)


# outlier detection
uhat <- (y - drop(X %*% fit$coef)) / fit$sigma
ns <- dim(sims)[1]
u <- matrix(nrow = ns, ncol = n)
yhat <- matrix(nrow = ns, ncol = n)
for (s in 1:ns) {
  yhat[s,] <- drop(X %*% sims[s, 1, 1:p]) # fitted values
  u[s,] <- (y - yhat[s,]) / sims[s, 1, p + 1] # standardized residuals
}

# Bayesian residual plots:
outlier_threshold <- function (ntests = 1, alpha = .05)
  -qnorm(.5 * (1 - (1 - alpha) ^ (1 / ntests)))
#k <- outlier_threshold(1) # very liberal
k <- outlier_threshold(n) # very conservative, close to Bonferroni

boxplot(u, outline = FALSE)
abline(h = 0, lty = 3); abline(h = c(-k, k), lty = 2, col = "red")

plot(c(yhat), c(u), pch = 20, col = "#dfdfdf1a",
     xlab = "Fitted values", ylab = "Standardized residuals")
points(drop(X %*% fit$coef), uhat, pch = 3)
abline(h = 0, lty = 3); abline(h = c(-k, k), lty = 2, col = "red")
lines(lowess(c(yhat), c(u)), col = "red")

prob_outlier <- apply(u, 2, function (ui) mean(abs(ui) > k))
names(prob_outlier) <- 1:n
idx <- order(prob_outlier, decreasing = TRUE)
print(prob_outlier[idx])

op <- par(mfrow = c(3, 3))
for (i in idx[1:9]) {
  plot_hist(u[, i], main = i)
  abline(v = uhat[i], lwd = 2)
  abline(v = c(-k, k), lty = 2, col = "red")
}
par(op)
