n <- 200
ts.1 <- data.frame(t=1:n, y=arima.sim(list(ar = (0.9)), n))
plot(ts.1$t, ts.1$y, type="l")
ccf(ts.1$y, ts.1$y, lag.max=100)
grid(col="grey")
ccf.out <- cross.correlate(ts.1, method="iccf", max.lag = 100)
lines(ccf.out$tau, ccf.out$ccf, type="l", col="red", lwd=4)


ts.2 <- data.frame(t=1:n, y=arima.sim(list(ar = (0.55)), n))
ts.2$y <- ts.2$y + filter(ts.1$y, dnorm(seq(-5,5,by=0.2)), 
                          sides=1, circular = TRUE)
ts.2$y <- ts.2$y + rnorm(n, sd=1)
plot(ts.1$t, ts.1$y, type="l")
lines(ts.2$t, ts.2$y, col="red")

result <- ccf(ts.1$y, ts.2$y, lag.max=100, plot=FALSE)
plot(result$lag, result$acf, type="h", bty="n", xlab="lag", ylab="CCF")

mask <- sample(1:n, 50, replace = FALSE)
mask <- sort(mask)
ts.1 <- ts.1[mask, ]

ccf.out <- cross.correlate(ts.1, ts.2, method="iccf", max.lag = 100, dtau=0.1)
grid(col="grey")
lines(ccf.out$tau, ccf.out$ccf, type="l", col="red", lwd=4)
