//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; # number of observations
  int<lower=0> K; # number of variables
  vector[N] y;
  
    // For more general linear regression
   matrix[N, K] X;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  sigma ~ cauchy(0,2.5);
  y ~ normal(alpha+ X * beta, sigma);
}

generated quantities{
  vector[N] y_rep; // create another y replicates of same sample size
  for(n in 1:N){
    y_rep[n] = normal_rng(alpha + X[n,] * beta, sigma);  // random generated or something, have to use a for-loop
  }
}
