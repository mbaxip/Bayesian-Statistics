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
  vector[N] y;
  
    // For more general linear regression
   vector[N] x;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real a;
  real b;
  real mu;
  real<lower=0> sigma;
  real<lower=0> tau;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  a ~ normal(0,10);
  b ~ normal(0,10);
  mu ~ normal(0,10);
  tau ~ cauchy(0,2.5);
  sigma ~ cauchy(0,2.5);
  //x ~ normal(mu, tau+sigma);
  x ~ normal(mu, tau+2*sigma);
  y ~ normal(a+b*mu, b*tau+sigma);
}

// generated quantities{
//   vector[N] y_rep; // create another 39 y replicates of same sample size
//   vector[N] x_rep;
//   for(n in 1:N){
//     y_rep[n] = normal_rng(a+b*mu, b*tau+sigma); // random generated or something, have to use a for-loop
//     x_rep[n] = normal_rng(mu, tau+sigma);
//   }
// }
