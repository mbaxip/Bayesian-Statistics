data {
  int<lower=0> J; // #experiments
  int<lower=0> n[J]; // #trials per experiment
  int<lower=0> y[J]; // #successes per experiment
}
parameters {
  vector<lower=0, upper=1>[J] theta; // prob of success per experiment
  real<lower=0, upper=1> rho; // rho = alpha / (alpha + beta)
  real<lower=1> kappa; // kappa = alpha + beta
}
model {
  y ~ binomial(n, theta);
  theta ~ beta(rho * kappa, (1 - rho) * kappa);
  // P(rho) propto 1
  kappa ~ pareto(1, 0.5);
}
