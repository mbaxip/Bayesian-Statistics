theta <- rgamma(n = 10000,shape = 6919,scale = 1/10)
pred_y <- rpois(n = 10000,lambda = theta)
sort_pred_y<- sort(pred_y)
conf_interval_95 <- sort_pred_y[c(251,9750)]
print(conf_interval_95)