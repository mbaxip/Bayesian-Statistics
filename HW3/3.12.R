#b
t = seq(1,10,by=1)
y = c(24,25,31,31,22,21,26,20,16,22)
mod <- lm(y~t)
summary(mod)
alpha_mean <- mod$coefficients[1]
alpha_sd <- coef(summary(mod))[1, 2]
beta_mean <- mod$coefficients[2]
beta_sd <- coef(summary(mod))[2, 2]
# alpha and beta follow independent normal distribution
alpha.grid=seq(alpha_mean-3*alpha_sd, alpha_mean+3*alpha_sd, length=100)
beta.grid=seq(beta_mean-3*beta_sd, beta_mean+3*beta_sd, length=100)

p=matrix(NA, 100, 100)
evaluate_dist=function(alpha, beta){
  product_alpha_beta_dist = dnorm(alpha, alpha_mean, alpha_sd)*dnorm(beta, beta_mean, beta_sd)
  return(product_alpha_beta_dist)
  }

for (i in 1:100){
  for (j in 1:100){
    p[i,j]=evaluate_dist(alpha.grid[i], beta.grid[j])}}
image(alpha.grid,beta.grid,p)
points(alpha_mean,beta_mean,pch=3)
contour(alpha.grid, beta.grid, p,xlab = "alpha",ylab = "beta",add=TRUE)

#f
#alpha.grid=seq(alpha_mean-3*alpha_sd, alpha_mean+3*alpha_sd, length=100)
#beta.grid=seq(beta_mean-3*beta_sd, beta_mean+3*beta_sd, length=100)

alpha.grid=seq(20, 40, length=100)
beta.grid=seq(-2.0, 0.5, length=100)
evaluate_postdist=function(alpha, beta){
  post_dist = exp(-(10*alpha+sum(beta*t)) + sum(y*log(alpha+beta*t)))
  return(post_dist)
}
z=matrix(NA, 100, 100)
for (i in 1:100){
  for (j in 1:100){
    z[i,j]=evaluate_postdist(alpha.grid[i], beta.grid[j])}}
image(alpha.grid,beta.grid,z,xlim=c(20, 40), ylim=c(-2.0, 0.5))
points(alpha_mean,beta_mean,pch=3)
contour(alpha.grid, beta.grid, z,xlab = "alpha",ylab = "beta",  add=TRUE)
#contour(alpha.grid, beta.grid, z, xlim=c(20, 40),
        #ylim=c(-2.5, 0.5), xlab="alpha", ylab="beta")

#1000 draws from posterior
pa <- rowSums(z) # marginal on alpha
ns <- 1000 # #samples
alpha.s <- beta.s <- numeric(ns)
for (s in 1:ns) {
  ia <- sample.int(100, 1, prob=pa) # sample alpha
  ib <- sample.int(100, 1, prob=z[ia,]) # sample beta | alpha
  alpha.s[s] <- alpha.grid[ia]; beta.s[s] <- beta.grid[ib]
}
# add to plot
points(alpha.s, beta.s, pch=19, col='blue', cex=.6, xlim=c(20, 40), ylim=c(-2.0, 0.5))


## g
t1 <- 11
post_theta <- alpha.s+ t1*beta.s
hist(post_theta, col = 'blue')

## h
postpred_y <- rpois(1000,alpha.s+ t1*beta.s)
sort_postpred <- sort(postpred_y)
conf_interval_95 <- sort_postpred[c(25,975)]
print(conf_interval_95)
# [1] 10 30



