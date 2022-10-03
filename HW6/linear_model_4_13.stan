// 14.13 stan file
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
  vector[N] x1;
  vector[N] x2;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real a;
  real b;
  real c;
  real mu1;
  real mu2;
  real<lower=0> sigma;
  real<lower=0> tau1;
  real<lower=0> tau2;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  a ~ normal(0,10);
  b ~ normal(0,10);
  c ~ normal(0,10);
  mu1 ~ normal(0,10);
  mu2 ~ normal(0,10);
  tau1 ~ cauchy(0,2.5);
  tau2 ~ cauchy(0,2.5);
  sigma ~ cauchy(0,2.5);
  x1 ~ normal(mu1, tau1+sigma);
  x2 ~ normal(mu2, tau2+sigma);
  y ~ normal(a+b*mu1+c*mu2, b*tau1+c*tau2+sigma);
}


