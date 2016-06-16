
# ----------------------------------------------
# load the data

ts.1 <- read.table("http://www.astronomy.ohio-state.edu/~agnwatch/n5548/lcv/c5100.dat")
ts.2 <- read.table("http://www.astronomy.ohio-state.edu/~agnwatch/n5548/lcv/hbeta.dat")
colnames(ts.1) <- c("t", "y", "dy")
colnames(ts.2) <- c("t", "y", "dy")

# plot the two time series
plot(ts.1[,1], ts.1[,2], type="o", bty="l")
segments(ts.1[,1], ts.1[,2]-ts.1[,3], ts.1[,1], ts.1[,2]+ts.1[,3])
lines(ts.2[,1], ts.2[,2], type="o", col="red")
segments(ts.2[,1], ts.2[,2]-ts.2[,3], ts.2[,1], ts.2[,2]+ts.2[,3], col="red")
grid(col="darkgrey")

# ----------------------------------------------
# compute the DCF
result <- cross.correlate(ts.1, ts.2, local.est=TRUE, zero.clip=TRUE,
                          dtau=1, nsim=500, max.lag=150)

# plot the DCF
plot(result$tau, result$ccf, type="o", bty="n", ylim=c(-1, 1),
     xlab="lag", ylab="DCF")
grid(col="darkgrey")

# ----------------------------------------------
# plot distributions of peaks and centroids

hist(result$cent.dist, breaks=40, col="steelblue1", prob=TRUE,
     border=NA, main="CCF centroid distribution")
lines(density(result$cent.dist, n=256, na.rm=TRUE), lwd=2)

# ----------------------------------------------
