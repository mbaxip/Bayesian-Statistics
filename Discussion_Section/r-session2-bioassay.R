# MA 578: R session, 09/24/19
# y_i | alpha, beta ~ Binomial(n_i, inv.logit(alpha + beta * x_i))
# P(alpha, beta) propto 1

# read data (Table 3.1 in BDA)
bioassay <- read.csv('data/bioassay.csv', header=TRUE)

# check frequentist estimates
gb <- glm(cbind(deaths, n - deaths) ~ logdose, family=binomial, data=bioassay)
summary(gb) # check posterior mode and standard deviation to define grid range

m <- 100 # number of grid subdivisions
alpha <- seq(-2, 5, length=m) # (-5, 10) in Gelman
beta <- seq(-5, 25, length=m) # (-10, 40) in Gelman

ilogit <- function (x) 1 / (1 + exp(-x))
lhood1 <- function (a, b) { # log-likelihood
  p <- ilogit(a + b * bioassay$logdose)
  sum(bioassay$deaths * log(p) + (bioassay$n - bioassay$deaths) * log(1 - p))
}

# (alpha, beta) ~ N(c(a.mean, b.mean), ab.prec^(-1))
a.mean <- 0; a.sd <- 2
b.mean <- 10; b.sd <- 10
ab.corr <- .5
ab.s2 <- ab.corr * a.sd * b.sd
ab.prec <- solve(matrix(c(a.sd^2, ab.s2, ab.s2, b.sd^2), nrow=2))
lprior <- function (a, b) {
  v <- matrix(c(a - a.mean, b - b.mean), ncol=1)
  return(-.5 * crossprod(v, ab.prec %*% v))
}

lj <- matrix(nrow=m, ncol=m) # log-joint grid
for (ia in 1:m)
  for (ib in 1:m)
    lj[ia, ib] <- lhood1(alpha[ia], beta[ib]) + lprior(alpha[ia], beta[ib])
pab <- exp(lj - log(sum(exp(lj)))) # posterior on alpha, beta
# plot
image(alpha, beta, pab)
points(coef(gb)[1], coef(gb)[2], pch=3)
contour(alpha, beta, pab, add=TRUE)

# sample alpha and beta
pa <- rowSums(pab) # marginal on alpha
ns <- 1000 # #samples
alpha.s <- beta.s <- numeric(ns)
for (s in 1:ns) {
  ia <- sample.int(m, 1, prob=pa) # sample alpha
  ib <- sample.int(m, 1, prob=pab[ia,]) # sample beta | alpha
  alpha.s[s] <- alpha[ia]; beta.s[s] <- beta[ib]
}
# add to plot
points(alpha.s, beta.s, pch=19, col='blue', cex=.6)

# plot distribution
quantile(alpha.s, c(.025, .975)) # 95% posterior interval
hist(alpha.s); abline(v=coef(gb)[1], lwd=2) # line marks mode
abline(v=quantile(alpha.s, c(.025, .975)), lty=2)
quantile(beta.s, c(.025, .975)) # 95% posterior interval
hist(beta.s); abline(v=coef(gb)[2], lwd=2)
abline(v=quantile(beta.s, c(.025, .975)), lty=2)


# LD50: log dose x such that logit(alpha + beta * x) = 0, i.e.,
# x = -alpha / beta
quantile(-alpha.s / beta.s, c(.025, .975)) # 95% posterior interval
hist(-alpha.s / beta.s, nclass=20); abline(v=-coef(gb)[1] / coef(gb)[2])
abline(v=quantile(-alpha.s / beta.s, c(.025, .975)), lty=2)


# [ EXTRA ]
# prediction at tilde(logdose)
x <- -.1
quantile(ilogit(alpha.s + beta.s * x), c(.025, .975))

logx <- seq(-1, 1, length=m)
li <- sapply(logx, function (x) quantile(ilogit(alpha.s + beta.s * x), .025))
mi <- sapply(logx, function (x) quantile(ilogit(alpha.s + beta.s * x), .5))
ui <- sapply(logx, function (x) quantile(ilogit(alpha.s + beta.s * x), .975))
plot(logx, mi, ylim=c(0,1), type='n',
     xlab='log dose', ylab=expression(P(tilde(y)~"="~1~"|"~y)))
for (i in 1:m) lines(c(logx[i], logx[i]), c(li[i], ui[i]))
points(logx, mi, pch=19, cex=.6)
ei <- sapply(logx, function (x) mean(ilogit(alpha.s + beta.s * x)))
points(logx, ei, pch=19, cex=.6, col='red')

