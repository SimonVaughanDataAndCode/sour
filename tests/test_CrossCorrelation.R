
# package the time series in data frames
source("cross_correlation.R")
ts.1 <- read.table("data/time_series_1.txt", header=TRUE)
ts.2 <- read.table("data/time_series_2.txt", header=TRUE)

# insert some dodgy (missing) data

ts.1[20:50,2] <- NA
ts.2[40:70,2] <- NA

plot(ts.1$t, ts.1$y, type="s", bty="n", xlab = "time", 
     ylab = "value")
lines(ts.2$t, ts.2$y, type="s", lwd=2, col="red")

# ----------------------------------------------
# compute the DCF
dcf.out <- cross.correlate(ts.1, ts.2, local.est = TRUE, 
                           dtau = 1, nsim = 2000, max.lag = 70)

dtau <- diff(dcf.out$tau[1:2])

# plot the CCF
plot(dcf.out$tau-dtau/2, dcf.out$ccf, type = "s", lwd = 2, 
     bty = "n", ylim = c(-1, 1), xlab = "lag", ylab = "CCF", 
     col = "black", main = "DCF test - true lag is 10")
grid(col = "darkgrey")

# plot distributions of centroids
hist(dcf.out$cent.dist, breaks = 50, col = "steelblue1", 
     border = NA, main = "CCF centroid distribution", prob = TRUE)
lines(density(dcf.out$cent.dist, n = 256, na.rm = TRUE), lwd = 2)
cat('-- mean lag', mean(dcf.out$cent.dist), fill = TRUE)
