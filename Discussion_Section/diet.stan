data {
  int<lower=0> J;
  int<lower=0> n[J];
  int<lower=0> idx[J]; // index shifts in y
  int<lower=0> N; // sum(n)
  real y[N];
}
parameters {
  vector[J] theta;
  real mu;
  real<lower=0> tau;
  real<lower=0> sigma;
}
model {
  for (j in 1:J)
    for (i in 1:n[j]) y[idx[j] + i] ~ normal(theta[j], sigma);
  for (j in 1:J) theta[j] ~ normal(mu, tau);
  target += -log(sigma);
}
