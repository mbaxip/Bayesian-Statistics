# MA 578: R session, 10/22/19
# Model checking

library(beanplot)

# [ Binomial hierarchical model (BDA Section 3.3) ]
# y_j | theta_j ~ind Binom(n_j, theta_j)
# theta_j | kappa ~iid Beta(r * kappa, (1 - r) * kappa)
# kappa ~ Gamma(beta, beta)

tumor <- read.csv("data/rattumor.csv", header = TRUE)
r <- with(tumor, sum(y) / sum(n)) # "empirical Bayes"
beta <- 10

# prior predictive
ns <- 1000
kappa_s <- rgamma(ns, beta, beta)
theta_s <- rbeta(ns, r * kappa_s, (1 - r) * kappa_s)
hist(theta_s, col = "gray", border = "white")
nm <- round(mean(tumor$n)) # "typical" experiment
y_s <- rbinom(ns, nm, theta_s) / nm
hist(y_s, col = "gray", border = "white")

lprior <- function (kappa) (beta - 1) * log(kappa) - beta * kappa
lhood_k <- function (k) { # log(P(y | kappa))
  a <- r * k; b <- (1 - r) * k
  with(tumor,
    sum(lgamma(a + y) - lgamma(a) + lgamma(b + n - y) - lgamma(b) -
        (lgamma(a + b + n) - lgamma(a + b))))
}

log_normalize <- function (x) {
  x <- x - max(x) # avoid underflows
  exp(x - log(sum(exp(x))))
}

m <- 100
k <- seq(0, 6, length = m + 1)[-1]
pp <- log_normalize(lprior(k))
pk <- log_normalize(sapply(k, function (x) lprior(x) + lhood_k(x)))
plot(k, pk, type = "l", ylim = c(0, max(pp, pk))) # posterior
lines(k, pp, lty = 2) # prior

# posterior predictive
ns <- 1000
k_s <- sample(k, ns, replace = TRUE, prob = pk)
theta_s <- matrix(nrow = ns, ncol = nrow(tumor))
yrep <- matrix(nrow = ns, ncol = nrow(tumor))
for (j in 1:nrow(tumor)) {
  theta_s[,j] <- with(tumor,
                      rbeta(ns, r * k_s + y[j], (1 - r) * k_s + n[j] - y[j]))
  yrep[,j] <- rbinom(ns, tumor$n[j], theta_s[,j])
}
yrep <- sweep(yrep, 2, tumor$n, `/`) # switch to proportions

p_data <- with(tumor, y / n)
boxplot(yrep, outline = FALSE, col = "gray")
points(p_data, pch = 19, cex = .8, col = "red")
abline(h = with(tumor, sum(y) / sum(n)), col = "red", lty = 2)

# overall
beanplot(list(y = tumor$y / tumor$n, y_pred = yrep), side = "both",
         col = list("black", "gray"), what = c(TRUE, TRUE, TRUE, FALSE))


# overall p-value:
# T(y, t) := |y^(.9) - t^(.5)| - |y^(.1) - t^(.5)|, t^(a) is a-quantile of t
tm <- apply(theta_s, 1, median) # t^(.5)
s <- quantile(p_data, c(.1, .9))
Tdata <- abs(s[2] - tm) - abs(s[1] - tm)
s <- abs(t(apply(yrep, 1, quantile, c(.1, .9))) - tm)
Trep <- s[,2] - s[,1]
mean(Trep >= Tdata)
beanplot(list(Tdata = Tdata, Trep = Trep), side = "both",
         col = list("black", "gray"), what = c(TRUE, TRUE, TRUE, FALSE))

# T(y, t) := sum|y - t|
Tdata <- abs(sweep(theta_s, 2, p_data, `-`))
Trep <- abs(yrep - theta_s)
mean(Trep >= Tdata)
beanplot(list(Tdata = Tdata, Trep = Trep), side = "both",
         col = list("black", "gray"), what = c(TRUE, TRUE, TRUE, FALSE))
