# MA 578: R session, 10/08/19
# Hierarchical Models

# [ Normal hierarchical model (BDA Section 5.5) ]
# y_.j | theta_j ~ind N(theta_j, sigma_j ^ 2)
# theta_j | mu, tau^2 ~iid N(mu, tau^2)

sat <- read.csv("data/sat.csv", header = TRUE)
y_mean <- with(sat, weighted.mean(mean, 1 / sd ^ 2))

zq <- qnorm(.75); delta <- .1 # plot params
mu0 <- 0; kappa0 <- 0 # prior params


mu_hat <- function (tau2)
  with(sat,
       (sum(mean / (tau2 + sd ^ 2)) + kappa0 * mu0 / tau2) /
       (sum(1 / (tau2 + sd ^ 2)) + kappa0 / tau2))

V_mu <- function (tau2)
  with(sat, 1 / (sum(1 / (tau2 + sd ^ 2)) + kappa0 / tau2))


log_Ptau <- function (tau) {
  tau2 <- tau ^ 2
  mu <- mu_hat(tau2)
  l <- with(sat, sum((mean - mu) ^ 2 / (sd ^ 2 + tau2))) +
          kappa0 * (mu0 - mu) ^ 2 / tau2
  lp <- ifelse(kappa0 > 0, -.5 * log(tau2), 0) + # log prior
              .5 * log(V_mu(tau2)) - .5 * sum(log(sat$sd ^ 2 + tau2)) - .5 * l
}

# plot P(tau | y)
m <- 100
tau <- seq(0, 30, length = m)[-1]
lt <- sapply(tau, log_Ptau)
ptau <- exp(lt - log(sum(exp(lt))))
plot(tau, ptau, type = "l", xlab = expression(tau),
     ylab = expression(P(tau~"|"~y)))
tau_m <- sum(tau * ptau) # point estimate: posterior mean

# First check: marginal samples
ns <- 1000 # #samples
tau2_s <- sample(tau, ns, replace = TRUE, prob = ptau) ^ 2 # tau^2 | y
mu_s <- rnorm(ns, sapply(tau2_s, mu_hat), sqrt(sapply(tau2_s, V_mu))) # mu | tau^2, y
thetat_s <- matrix(nrow = ns, ncol = nrow(sat))
ytilde_s <- matrix(nrow = ns, ncol = nrow(sat))
for (j in 1:nrow(sat)) {
  vj <- 1 / (1 / sat$sd[j] ^ 2 + 1 / tau2_s)
  tj <- (sat$mean[j] / sat$sd[j] ^ 2 + mu_s / tau2_s) * vj
  thetat_s[, j] <- rnorm(ns, tj, sqrt(vj)) # theta_j | mu, tau^2, y
  ytilde_s[, j] <- rnorm(ns, thetat_s[,j], sat$sd[j]) # y~ | theta_j, mu, tau^2, y
}

hist(mu_s, col = "gray", border = "white", xlab = expression(mu),
     prob = TRUE, main = expression(mu~"|"~y))
abline(v = mean(mu_s), lwd = 2) # MC-estimated posterior mean
abline(v = y_mean, lwd = 2, col = "red")

boxplot(thetat_s, names = sat$school, col = "gray", outline = FALSE,
        main = expression(theta[j]~"|"~y))
points(sat$mean, pch = 19, col = "red")
abline(h = y_mean, lty = 2, lwd = 2, col = 2)


# [ conditional posterior theta_j | tau2, y ]
theta_s <- function (tau2)
  with(sat, (mean / sd ^ 2 + mu_hat(tau2) / tau2) / (1 / sd ^ 2 + 1 / tau2))

V_s <- function (tau2) {
  V <- with(sat, 1 / (1 / sd ^ 2 + 1 / tau2))
  V + (V / tau2) ^ 2 * V_mu(tau2)
}

# [ plots for E[theta | tau, y] and sd[theta | tau, y], as in BDA
t <- matrix(nrow = length(tau), ncol = nrow(sat))
for (i in 1:length(tau)) t[i,] <- theta_s(tau[i] ^ 2)
plot(tau, rowMeans(t), type = "n", ylim = range(sat$mean),
     xlab = expression(tau), ylab = expression(E(theta[j]~"|"~tau, y)))
for (i in 1:nrow(sat)) lines(tau, t[,i], col = i)
abline(v = tau_m, lty = 2) # posterior mean for tau

for (i in 1:length(tau)) t[i,] <- sqrt(V_s(tau[i] ^ 2))
plot(tau, rowMeans(t), type = "n", ylim = range(t),
     xlab = expression(tau), ylab = expression(sd(theta[j]~"|"~tau, y)))
for (i in 1:nrow(sat)) lines(tau, t[,i], col = i)
abline(v = tau_m, lty = 2) # posterior mean for tau


# [ posterior predictive, marginal and conditional on mu ]
tau_s <- 20 # large enough to be different from tau_m
mu <- seq(min(ytilde_s), max(ytilde_s), length = m)
mu_m <- sapply(tau ^ 2, mu_hat)
mu_sd <- sqrt(sapply(tau ^ 2, V_mu))
mu_p <- sapply(mu, function (x) sum(dnorm(x, mu_m, mu_sd) * ptau))
mu_cp <- dnorm(mu, mu_hat(tau_s ^ 2), sqrt(V_mu(tau_s ^ 2)))
# "pooled" y~, compare to posterior mu and cond posterior:
hist(ytilde_s, prob = TRUE, col = "gray", border = "white",
     xlab = expression(tilde(y)), main = expression(tilde(y)~"|"~y),
     ylim = c(0, max(mu_p, mu_cp)))
lines(mu, mu_p, lwd = 2) # marginal posterior
lines(mu, mu_cp, lwd = 2, lty = 2) # cond posterior
abline(v = y_mean, lty = 2, col = 2)


# [ posterior on theta, theta | tau, y, compare to cond posterior with larger tau^2 ]
boxplot(thetat_s, names = sat$school, outline = FALSE, col = "gray", # marginal
        ylim = c(min(sat$mean - sat$sd), max(sat$mean + sat$sd)),
        main = expression(theta[j]~"|"~y))
t <- theta_s(tau_s ^ 2) # cond posterior
s <- sqrt(V_s(tau_s ^ 2))
points(1:nrow(sat) - delta, t, pch = 19)
points(1:nrow(sat) + delta, sat$mean, pch = 19, col = "red")
for (i in 1:nrow(sat)) {
  lines(c(i - delta, i - delta),
        c(t[i] - zq * s[i], t[i] + zq * s[i]), lwd = 3)
  lines(c(i + delta , i + delta),
        with(sat, c(mean[i] - zq * sd[i], mean[i] + zq * sd[i])),
        lwd = 3, col = "red")
}
abline(h = y_mean, lty = 2, lwd = 2, col = "red")


# [ posterior predictive: y~ | tau, y ~ N(theta_s, sd^2 + V_s*) ]
boxplot(ytilde_s, names = sat$school, outline = FALSE, col = "gray", # marginal
        main = expression(tilde(y)[j]~"|"~y))
t <- theta_s(tau_s ^ 2)
s <- sqrt(sat$sd ^ 2 + V_s(tau_s ^ 2))
points(1:nrow(sat) - delta, t, pch = 19)
points(1:nrow(sat) + delta, sat$mean, pch = 19, col = 2)
for (i in 1:nrow(sat)) {
  lines(c(i - delta , i - delta),
        c(t[i] - zq * s[i], t[i] + zq * s[i]), lwd = 3)
  lines(c(i + delta , i + delta),
        with(sat, c(mean[i] - zq * sd[i], mean[i] + zq * sd[i])),
        lwd = 3, col = 2)
}
abline(h = y_mean, lty = 2, lwd = 2, col = 2)



# [ Binomial hierarchical model (BDA Section 5.3) ]
# y_j | theta_j ~ind Binom(n_j, theta_j)
# theta_j | alpha, beta ~iid Beta(alpha, beta)

tumor <- read.csv("data/rattumor.csv", header = TRUE)

lprior <- function (a, b) -2.5 * log(a + b)
lhood_ab <- function (a, b)
  with(tumor,
    sum(lgamma(a + y) - lgamma(a) + lgamma(b + n - y) - lgamma(b) -
        (lgamma(a + b + n) - lgamma(a + b))))

# [ numerical approximation to P(alpha, beta | y) ]
# find 'grid' boundaries; first find approximate MLE (am, bm):
# with lgamma(a + y) - lgamma(a) ~= y * log(a) we have
# am ~= r * ks and bm ~= (1 - r) * ks, with
r <- with(tumor, sum(y) / sum(n))
# 'ks' can be found numerically:
m <- 100
k <- seq(10, 20, length = m)
ks <- k[which.max(sapply(k, function (x) lhood_ab(r * x, (1 - r) * x)))]
# check maximum:
plot(k, sapply(k, function (x) lhood_ab(r * x, (1 - r) * x)), type = "l",
     ylab = "log-likelihood")
abline(v = ks, lty = 2)
# thus:
am <- r * ks; bm <- (1 - r) * ks

a <- seq(0, 2 * am, length = m + 1)[-1]
b <- seq(0, 2 * bm, length = m + 1)[-1]

lab <- matrix(nrow = m, ncol = m)
for (ia in 1:m)
  for (ib in 1:m)
    lab[ia, ib] <- lprior(a[ia], b[ib]) + lhood_ab(a[ia], b[ib])
pab <- exp(lab - log(sum(exp(lab))))
# posterior mode
im <- which(pab == max(pab), arr = TRUE)
as <- a[im[1]]; bs <- b[im[2]]
image(a, b, pab); points(as, bs, pch = 3)
contour(a, b, pab, add = TRUE)


# [BDA] with transformation: r(a, b) = log(a / b) and s(a, b) = log(a + b)
r <- seq(-2.2, -1.4, length = m) # as in BDA
s <- seq(1, 5, length = m) # as in BDA

lrs <- matrix(nrow = m, ncol = m)
for (ir in 1:m) {
  for (is in 1:m) {
    ilogit_r <- 1 / (1 + exp(-r[ir]))
    exp_s <- exp(s[is])
    a_rs <- exp_s * ilogit_r
    b_rs <- exp_s * (1 - ilogit_r)
    log_jacobian <- 2 * s[is] + log(ilogit_r * (1 - ilogit_r))
    lrs[ir, is] <- log_jacobian + lprior(a_rs, b_rs) + lhood_ab(a_rs, b_rs)
  }
}
prs <- exp(lrs - log(sum(exp(lrs))))
# posterior mode
im <- which(prs == max(prs), arr = TRUE)
rs <- r[im[1]]; ss <- s[im[2]]
image(r, s, prs); points(rs, ss, pch=3)
contour(r, s, prs, add = TRUE)


# sampling alpha, beta, and theta
pa <- rowSums(pab)
ind_ab <- matrix(1:(m ^ 2), nrow = m)
sample_ab <- function (ns) {
  is <- sapply(sample.int(m ^ 2, ns, replace = TRUE, prob = pab),
               function (ind) which(ind_ab == ind, arr = TRUE))
  cbind(a[is[1,]], b[is[2,]])
}

ns <- 1000
ab_s <- sample_ab(ns)
theta_s <- matrix(nrow = ns, ncol = nrow(tumor))
ytilde_s <- matrix(nrow = ns, ncol = nrow(tumor))
for (j in 1:nrow(tumor)) {
  theta_s[,j] <- with(tumor, rbeta(ns, ab_s[,1] + y[j], ab_s[,2] + n[j] - y[j]))
  ytilde_s[,j] <- with(tumor, rbinom(ns, n[j], theta_s[,j]))
}


# [ conditional posterior on theta ]
# numerical approximation based on alpha,beta-grid
t <- seq(0, 1, length = m + 1)[-1] # note: avoid 0
theta_p <- matrix(nrow = m, ncol = nrow(tumor))
for (j in 1:nrow(tumor)) {
  for (i in 1:m) {
    theta_p[i, j] <- sum(pab * outer(a, b,
        function (a, b) with(tumor, dbeta(t[i], a + y[j], b + n[j] - y[j]))))
  }
}

# plot each theta_j marginally
j <- sample.int(nrow(tumor), 1)
tlj <- with(tumor, dbeta(t, 1 + y[j], 1 + n[j] - y[j])) # likelihood
tcpj <- with(tumor, dbeta(t, as + y[j], bs + n[j] - y[j])) # cond posterior
# predictive proportions:
hist(ytilde_s[,j] / tumor$n[j],
     border = "white", col = "gray", prob = TRUE,
     xlim = c(0, 1), ylim = c(0, max(tlj, tcpj, theta_p[,j])),
     xlab = expression(theta[j]), ylab = "density",
     main = paste("[ j = ", j, "] y / n = ", tumor$y[j], "/", tumor$n[j]))
lines(t, tlj, col = "blue") # likelihood
lines(t, tcpj, lty = 2, lwd = 2) # cond posterior
lines(t, theta_p[,j], lwd = 2) # marginal posterior
abline(v = with(tumor, y[j] / n[j]), col = "blue")
abline(v = with(tumor, sum(y) / sum(n)), col = "red", lty = 2)


theta_l <- apply(theta_s, 2, quantile, .025)
theta_q1 <- apply(theta_s, 2, quantile, .25)
theta_med <- apply(theta_s, 2, quantile, .5)
theta_q3 <- apply(theta_s, 2, quantile, .75)
theta_u <- apply(theta_s, 2, quantile, .975)

at <- seq(1, length = nrow(tumor), by = 1.5)
boxplot(sweep(ytilde_s, 2, tumor$n, `/`),
        outline = FALSE, col = "gray", at = at,
        main = "Posterior predictive", ylab = "proportion")
points(at, tumor$y / tumor$n, pch = 19, cex = .8, col = "red") # data
delta <- .4
points(at + delta, theta_med, pch = 19, cex = .8, col = "blue")
for (j in 1:nrow(tumor)) { # marginal posterior
  lines(c(at[j] + delta, at[j] + delta), c(theta_q1[j], theta_q3[j]),
        lwd = 3, col = "blue")
  lines(c(at[j] + delta, at[j] + delta), c(theta_l[j], theta_u[j]),
        col = "blue")
}
abline(h = with(tumor, sum(y) / sum(n)),
       col = "red", lwd = 2, lty = 2)

