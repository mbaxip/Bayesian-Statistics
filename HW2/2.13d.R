theta <- rgamma(n=10000,shape = 6919,scale = 1/(5.715*10^12))
y_pred <- rpois(n=10000,lambda = 8*10^11 * theta)
sort_pred_y <- sort(y_pred)
conf_interval_95 <- sort_pred_y[c(251,9750)]
print(conf_interval_95)