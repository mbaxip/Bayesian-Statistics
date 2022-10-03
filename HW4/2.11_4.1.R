y = c(-2, -1, 0, 1.5, 2.5)
# θ uniform on [−4, 4].
theta = seq(-4,4,by=0.001) #theta grid

# Unormalized Posterior function
unnorm_posterior <- function (theta) {
  unnorm_post = prod(1/(1+(y-theta)^2))
  return(unnorm_post)
}
m = length(theta)
# Unnormalized Posterior grid
unnorm_poster <- matrix(0,m,1)
for (i in 1:m)
  unnorm_poster[i] <- unnorm_posterior(theta[i])

plot(theta, unnorm_poster, title('UnNormalized Posterior Density Function'), xlab = 'Theta', ylab = 'p(Theta|Y)', col = 'red')

# Normalized Posterior grid
prob_y = integrate(Vectorize(unnorm_posterior),-Inf,Inf)
norm_post <- unnorm_poster/prob_y$value

# 2.11 a
plot(theta, norm_post, title('Normalized Posterior Density Function'), xlab = 'Theta', ylab = 'p(Theta|Y)', col = 'red')

# 4.1 b

# First Derivative of log Posterior function
FirstDeriv_lPost <- function(theta) {
  first_Deriv <- 2 * sum((y-theta)/(1+(y-theta)^2))
  return(first_Deriv)
}

# Iteratively solving for Posterior Mode

NewtonRhapson <- uniroot(FirstDeriv_lPost, c(-4,4))

Posterior_Mode <- NewtonRhapson$root #[1] 0.08759136

## 4.1 c
# Fisher Information (Posterior Mode) : I_PostMode
I_PostMode_inv = 1/(2 * sum ((1-(y-Posterior_Mode)^2)/(1+(y-Posterior_Mode)^2)^2))
Norm_Approx <- dnorm(theta, mean = Posterior_Mode, sd = sqrt(I_PostMode_inv))
plot(theta, Norm_Approx, title('Normal Approximation of Unnormalized Posterior Density'), xlab = 'Theta', ylab = 'Norm Approx: p(Theta|Y)', col = 'blue')


  
  
  



