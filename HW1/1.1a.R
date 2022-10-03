y <- seq(-10,10,.001) 
prob_y <- 0.5*(dnorm(y,1,2) + dnorm(y,2,2)) # pdf of y formula
plot (y, prob_y,type='l')