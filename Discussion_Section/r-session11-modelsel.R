# MA 578 - Bayesian Statistics
# Model selection in linear models with semi-conjugate prior
# y | beta ~ N(X * beta, sigma2 * I_n)
# beta | sigma2 ~ N(beta0, sigma2 * Sigma0)
# sigma2 ~ Inv-Chisq(nu, tau2)

# y ~ t(nu, X * beta0, tau2 * (I_n + X * Sigma0 * X'))
log_marg_y <- function (y, x, beta0, Sigma0, nu, nu_tau2) {
  n <- length(y)
  C <- chol(diag(n) + x %*% tcrossprod(Sigma0, x))
  v <- crossprod(backsolve(C, y - x %*% as.matrix(beta0), transpose = TRUE))
  -sum(log(diag(C))) - (nu + n) / 2 * log(nu_tau2 + v)
}


# Table 7.3 in BDA: radon measurements for three counties in MN
radon <- read.csv("data/radon.csv", comment = "#")
y <- radon$radon
x <- model.matrix(~ county + firstfloor + county:firstfloor, data = radon)
p <- ncol(x); n <- length(y)
nu <- tau2 <- 0 # non-informative for sigma2
beta0 <- rep(0, p); Sigma0 <- chol2inv(chol(crossprod(x))) # Zellner's g-prior

# test for no interaction
x0 <- model.matrix(~ county + firstfloor, data = radon)
p0 <- ncol(x0); ind0 <- 1:p0
Sigma0_null <- chol2inv(chol(crossprod(x0)))

# compute Bayes factors for a range of 'g'
m <- 100
g <- seq(0, n, length=m + 1)[-1]
BF <- numeric(m)
for (i in 1:m) {
  l1 <- log_marg_y(y, x, beta0, g[i] * Sigma0, nu, nu * tau2) # log(P(y|z_1))
  l0 <- log_marg_y(y, x[, ind0], beta0[ind0], g[i] * Sigma0_null,
                   nu, nu * tau2) # log(P(y|z_0))
  BF[i] <- exp(l1 - l0)
}

plot(g, BF, type = "l")
gamma <- 1; alpha <- .5
pH0 <- 1 - alpha
abline(h = gamma * pH0 / (1 - pH0), col = "red")
abline(v = g[which.max(bf)], lty = 2)


# [ multiple tests ]
#      1 county firstfloor county:firstfloor
z <- c(1, 0, 0,          0,            0, 0,
       1, 1, 1,          0,            0, 0,
       1, 0, 0,          1,            0, 0,
       1, 1, 1,          1,            0, 0,
       1, 1, 1,          1,            1, 1)
z <- matrix(z, ncol = ncol(x), byrow = TRUE)

g <- n; gamma <- 1; alpha <- .5
logit_a <- log(alpha / (1 - alpha))

pz <- apply(z, 1, function (zz) {
            ind_z <- which(zz == 1)
            Sigma0_z <- drop(g * chol2inv(chol(crossprod(x[, ind_z]))))
            log_marg_y(y, x[, ind_z], beta0[ind_z], Sigma0_z,
                       nu, nu * tau2) + sum(zz) * logit_a
       })
pz <- pz - max(pz); pz <- exp(pz - log(sum(exp(pz)))) # P(z | y)
pz1 <- colSums(sweep(z, 1, pz, `*`)) # P(z_j = 1 | y)



# http://www.stat.washington.edu/people/pdhoff/Book/Data/data/yX.o2uptake
o2uptake <- structure(c(-0.87, -10.74, -3.27, -1.97, 7.5, -7.25, 17.05, 4.96, 
10.4, 11.05, 0.26, 2.51, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 23, 22, 22, 25, 27, 20, 31, 
23, 27, 28, 22, 24, 0, 0, 0, 0, 0, 0, 31, 23, 27, 28, 22, 24), .Dim = c(12L, 
5L), .Dimnames = list(NULL, c("uptake", "intercept", "aerobic", 
"age", "aerobic.age")))
y <- o2uptake[, 1]; x <- o2uptake[, -1]

#      1 aerobic age aerobic:age 
z <- c(1,     0,  0,          0,
       1,     1,  0,          0,
       1,     0,  1,          0,
       1,     1,  1,          0,
       1,     1,  1,          1)
z <- matrix(z, ncol = ncol(x), byrow = TRUE)

p <- ncol(x); n <- length(y)
nu <- tau2 <- 0 # non-informative for sigma2
beta0 <- rep(0, p)

g <- n; gamma <- 1; alpha <- .5
logit_a <- log(alpha / (1 - alpha))

pz <- apply(z, 1, function (zz) {
            ind_z <- which(zz == 1)
            Sigma0_z <- drop(g * chol2inv(chol(crossprod(x[, ind_z]))))
            log_marg_y(y, x[, ind_z], beta0[ind_z], Sigma0_z,
                       nu, nu * tau2) + sum(zz) * logit_a
       })
pz <- pz - max(pz); pz <- exp(pz - log(sum(exp(pz)))) # P(z | y)
pz1 <- colSums(sweep(z, 1, pz, `*`)) # P(z_j = 1 | y)
