n <- 400
setwd("C:/Users/sav2/Documents/R/sour/tests")

# ---------------------------------------------------------
# white noise 1
ts.1 <- data.frame(t=1:n, y=arima.sim(list(ar = (0.0)), n))

# white noise 2
ts.2 <- data.frame(t=1:n, y=arima.sim(list(ar = (0.0)), n))

write.table(ts.1, "w1.txt", row.names = FALSE)
write.table(ts.2, "w2.txt", row.names = FALSE)

# random subset 1
mask <- sample(1:n, 100, replace = FALSE)
mask <- sort(mask)
ts.1r <- ts.1[mask, ]

# random subset 2
mask <- sample(1:n, 50, replace = FALSE)
mask <- sort(mask)
ts.2r <- ts.2[mask, ]

write.table(ts.1r, "w1r.txt", row.names = FALSE)
write.table(ts.2r, "w2r.txt", row.names = FALSE)
# ---------------------------------------------------------
# red  noise 1
ts.1 <- data.frame(t=1:n, y=arima.sim(list(ar = (0.95)), n))
ts.1n <- ts.1
ts.1n$y <- ts.1n$y + rnorm(n)

# white noise 2
ts.2 <- data.frame(t=1:n, y=arima.sim(list(ar = (0.95)), n))

write.table(ts.1, "r1.txt", row.names = FALSE)
write.table(ts.1n, "r1n.txt", row.names = FALSE)
write.table(ts.2, "r2.txt", row.names = FALSE)

# random subset 1
mask <- sample(1:n, 100, replace = FALSE)
mask <- sort(mask)
ts.1r <- ts.1[mask, ]

# random subset 1n
ts.1rn <- ts.1n[mask, ]

# random subset 2
mask <- sample(1:n, 50, replace = FALSE)
mask <- sort(mask)
ts.2r <- ts.2[mask, ]

write.table(ts.1r, "r1r.txt", row.names = FALSE)
write.table(ts.1rn, "r1nr.txt", row.names = FALSE)
write.table(ts.2r, "r2r.txt", row.names = FALSE)

# ---------------------------------
# convolution

ts.2 <- data.frame(t=1:n, y=arima.sim(list(ar = (0.85)), n))
ts.2$y <- ts.2$y + filter(ts.1$y, dnorm(seq(-5,5,by=0.2)), 
                          sides=1, circular = TRUE)
ts.2$y <- ts.2$y + rnorm(n, sd=1)
write.table(ts.2, "r2c.txt", row.names = FALSE)

ts.2r <- ts.2[mask, ]
write.table(ts.2r, "r2cr.txt", row.names = FALSE)
