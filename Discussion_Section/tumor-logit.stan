data {
  int<lower=0> J; // #experiments
  int<lower=0> n[J]; // #trials per experiment
  int<lower=0> y[J]; // #successes per experiment
}
parameters {
  vector<lower=0, upper=1>[J] theta; // prob of success per experiment
  real rho; // rho = logit(alpha / (alpha + beta))
  real kappa; // kappa = log(alpha + beta)
}
model {
  y ~ binomial(n, theta);
  theta ~ beta(inv_logit(rho) * exp(kappa), (1 - inv_logit(rho)) * exp(kappa));
  rho ~ logistic(0, 1);
  kappa ~ exponential(.5);
}
