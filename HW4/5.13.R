library(rootSolve)

y <- c(16, 9, 10, 13, 19, 20, 18, 17, 35, 55)
n <- c(74, 99, 58, 70, 122, 77, 104, 129, 308, 119)

## Estimate theta from data (Bin distribution)
theta_hat <- y/n
estimate_mean <- mean(theta_hat) #0.1961412
estimate_var <- var(theta_hat) #0.01112067

## Estimate alpha and beta from estimated mean and var for theta (Beta Distribution)
#par[1] is alpha and par[2] is beta
alpha_beta_Estimate          <- function(par) {
  diff  <- numeric(2)
  diff[1] <- par[1]/(par[1]+par[2]) - mean(y/n)
  diff[2] <- par[1]*par[2]/(((par[1]+par[2])^2)*(par[1]+par[2]+1)) - sd(y/n)^2
  diff
}
Estimate <- multiroot(alpha_beta_Estimate,c(1,1)) 
alpha_estimate <- round(Estimate$root[1],1)
beta_estimate <- round(Estimate$root[2],1)
#res1 <- paste("(",alpha_estimate,",",beta_estimate,")",sep=""):  "(2.6,10.6)"

## Creating a alpha-beta grid around the estimates from above

marginal_poster <-   function(alpha, beta) {
  post <-  1
  for (i in 1:length(y)) {
    if (n[i] > 100) n[i] = 100
    post  = post * ( ( ( gamma(alpha + beta) ) / ( gamma(alpha) * gamma(beta) ) ) * 
                       ( ( gamma(alpha + y[i] ) * gamma(beta + n[i] - y[i]) ) / ( gamma(alpha + beta + n[i]) ) ) )
  }
  # The hyper prior is defined below
  Hyper_prior(alpha,beta) * post
}

Hyper_prior <-  function(alpha,beta) {
  alpha * beta * (alpha + beta)^(-5/2)
}

axis_x  <-  seq(log(alpha_estimate/beta_estimate)*1.5,
                        log(alpha_estimate/beta_estimate)/1.5,length.out =151) # log(alpha/beta)
axis_y <-  seq(log(alpha_estimate+beta_estimate)/2.5,
                        log(alpha_estimate+beta_estimate)*1.5,length.out =151) # log(alpha + beta)
beta            <-  exp(axis_y)/(exp(axis_x)+1)
alpha           <-  exp(axis_y+axis_x)/(exp(axis_x)+1)

log_marginal <- function(x1,x2){
  log(marginal_poster(x1, x2))
}
posterior_dens       <-  outer(alpha,beta,log_marginal)
posterior_dens       <-  exp(posterior_dens - max(posterior_dens))
posterior_dens       <-  posterior_dens/sum(posterior_dens)

contours        <-  seq(min(posterior_dens), max(posterior_dens), length=10)
contour(axis_x, axis_y, posterior_dens,levels=contours, xlab=expression( log(alpha/beta) ), 
        ylab=expression( log(alpha+beta) ), xlim=c( min( axis_x  ), max( axis_x  ) ) , 
        ylim=c( min( axis_y ), max( axis_y ) ), 
        drawlabels=FALSE, main="Joint posterior density: p(alpha,beta|y)")

## Draw Samples from p(alpha,beta|y)

samples <- 1000
# Sum over all beta to get the marginal of alpha
marginal_alpha_dens  <-  apply(posterior_dens ,1, sum)
# sample_log_x: log(alpha/beta)
sample_log_x  <-  sample(axis_x, samples, replace=TRUE, prob = marginal_alpha_dens) 

# Compute conditional probability (p(beta|alpha))
conditional_prob_beta  <-  function(x) 
{
  posterior_dens[which(axis_x == sample_log_x[x]),]
}
# Sample beta according the the conditional probatility above
#sample_log_y: log(alpha + beta)
sample_log_y <-  sapply(1:samples,function(x) sample(axis_y,1,replace=TRUE,prob=conditional_prob_beta (x)))

# Add a uniform random jitter centered at zero with width equal to the grid spacing to make 
# simulation draws more continuous.   
grid.alpha         <-  axis_x[2] - axis_x[1]
grid.beta         <-  axis_y[2] - axis_y[1]
sample_log_y            <-  sample_log_y + runif(length(sample_log_y),-grid.beta/2,grid.beta/2)
sample_log_x            <-  sample_log_x + runif(length(sample_log_x),-grid.alpha/2,grid.alpha/2)

# Plot the sampled values
points(sample_log_x, sample_log_y,col = 'red', xlab=expression( log(alpha/beta)^s ), 
       ylab=expression( log(alpha+beta)^s ), xlim=c( min(axis_x) , max(axis_x) ) , 
       ylim=c( min(axis_y) , max(axis_y) ), 
       main="Sample Draws of log(alpha/beta) and log(alpha+beta)")

sample_beta          <-  exp( sample_log_y ) / ( exp(sample_log_x)+1 )
sample_alpha         <-  exp( sample_log_y + sample_log_x ) / ( exp(sample_log_x)+1 ) 


## 5.13c
# For each draw of hyper-parameters, draw a sample of θ from p(θ|alpha,beta,y)

theta_dist  <-  sapply(1:10,  function(x) rbeta(1000, sample_alpha +y[x], sample_beta + n[x] - y[x]))
theta_dist      <-  apply(theta_dist,2,sort)
plot(0:600/1000, 0:600/1000,  type="l", xlab="Observed rate",
     ylab="95% CI and median of posterior")

jitter.x        <-  y/n + runif(length(y),-0.01,0.01)
points(jitter.x, theta_dist[500,])
segments(jitter.x,theta_dist[25,], jitter.x,theta_dist[975,] )
title(main="Posterior Distribution of Bike rates for all 10 streets")


## 5.13d
# 1000 draws from Beta(sample_alpha,sample_beta):
  
sample_theta <- rbeta(1000, shape1 =sample_alpha , shape2 = sample_beta)   
CI <- round(sample_theta[order(sample_theta)][c(25,975)],2)

# The posterior interval for θ̂= (0.02, 0.46)
