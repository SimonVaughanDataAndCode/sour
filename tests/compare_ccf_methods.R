# ----------------------------------
# simulate two example time series
# The first is an AR1 process, the second is a filtered version of the first.

n <- 250
t <- 1:n
x <- arima.sim(model=list(ar=c(0.9)), n)
#y <- filter(x, filter=dchisq(seq(0,21,length=22), df=2), side=1, circular=TRUE)
y <- filter(x, filter=dnorm(seq(-10,10,length=21),mean=0.2143, sd=3.4), side=1, circular=TRUE)
y <- filter(y, filter=dnorm(seq(-10,10,length=21), sd=1.4), side=2, circular=TRUE)
y <- y + rnorm(n)*0.5
ts.1 <- data.frame(t=t, y=as.vector(x))
ts.2 <- data.frame(t=t, y=as.vector(y))

plot(ts.1, type="o", bty="n", pch=16)
points(ts.2, col="red", type="o")

# ----------------------------------
cat('-- DCF on good data.', fill=TRUE)

dcf.out <- cross.correlate(ts.1, ts.2, dtau=1, max.lag=50, nsim=5000, 
                           chatter=0)
dtau <- diff(dcf.out$tau[1:2])

plot(dcf.out$tau-dtau/2, dcf.out$ccf, type="s", lwd=2, bty="n", ylim=c(-1, 1),
     xlab="lag", ylab="CCF", col="black", main="CCF comparison\nTrue lag is -10.0")
grid(col="darkgrey")

# plot 90% band from sims
#pnk <- rgb(255, 192, 203, 100, maxColorValue = 255)
#polygon(c(dcf.out$tau, rev(dcf.out$tau)), c(dcf.out$lower, rev(dcf.out$upper)), 
#        col=pnk, border=NA)
#lines(dcf.out$tau - dtau/2, dcf.out$ccf, lwd=3, type="s")

# plot distributions of centroids
hist(dcf.out$cent.dist, breaks=50, col="steelblue1", 
     border=NA, main="CCF centroid distribution", prob=TRUE)
lines(density(dcf.out$cent.dist, n=256, na.rm=TRUE), lwd=2)
cat('-- mean lag', mean(dcf.out$cent.dist), fill=TRUE)

stop()

# plot 1 sigma band from Bartlett formula
dcf.out <- cross.correlate(ts.1, ts.2, dtau=1, max.lag=50, method="dcf", nsim=0,
                           chatter=0)
lines(dcf.out$tau, dcf.out$lower, lty=3, lwd=2)
lines(dcf.out$tau, dcf.out$upper, lty=3, lwd=2)

cat('-- ICCF on good data.', fill=TRUE)
iccf.out <- cross.correlate(ts.1, ts.2, max.lag=20, local.est=FALSE, dtau=0.25,
                            one.way=FALSE, method="iccf", nsim=500, chatter=0)
lines(iccf.out$tau, iccf.out$ccf, col="red", lwd=2)
# plot 'errors' too
lines(iccf.out$tau, iccf.out$lower, lty=3, lwd=2, col="red")
lines(iccf.out$tau, iccf.out$upper, lty=3, lwd=2, col="red")

ccf.out <- ccf(x, y, plot=FALSE, lag.max=100)
lines(ccf.out$lag, ccf.out$acf, col="blue", lwd=2)

legend(-50, -0.5, 
       c("CCF","DCF", "ICCF"),
       lty=c(1,1,1),
       lwd=c(2.5,2.5,2.5),col=c("blue","black","red")) 


# plot distributions of centroids
hist(iccf.out$cent.dist, breaks=50, col="steelblue1", 
     border=NA, main="CCF centroid distribution", prob=TRUE)
lines(density(iccf.out$cent.dist, n=256, na.rm=TRUE), lwd=2)

stop()

# ----------------------------------
# test ACFs

cat('-- DCF/ACF on good data.', fill=TRUE)
dcf.out <- cross.correlate(ts.1, dtau=1, max.lag=50, nsim=100)
dtau <- diff(dcf.out$tau[1:2])
plot(dcf.out$tau-dtau/2, dcf.out$ccf, type="s", lwd=2, bty="n", ylim=c(-1, 1),
     xlab="lag", ylab="CCF", col="black", main="ACF comparison")
grid(col="darkgrey")

# now swap ACF and compare
lines(-(dcf.out$tau-dtau/2), dcf.out$ccf, col="pink", lty=3, type="s")
lines(dcf.out$tau, dcf.out$lower, lty=3, lwd=2)
lines(dcf.out$tau, dcf.out$upper, lty=3, lwd=2)

cat('-- ICCF/ACF on good data.', fill=TRUE)
iccf.out <- cross.correlate(ts.1, method="iccf", max.lag=50, nsim=100)
lines(iccf.out$tau, iccf.out$ccf, col="pink", lwd=3)
lines(-iccf.out$tau, iccf.out$ccf, col="black", lwd=2, lty=3)
lines(iccf.out$tau, iccf.out$lower, lty=3, lwd=2, col="red")
lines(iccf.out$tau, iccf.out$upper, lty=3, lwd=2, col="red")

ccf.out <- ccf(x, x, plot=FALSE, lag.max=100)
lines(ccf.out$lag, ccf.out$acf, col="blue", lwd=2)

legend(-50, -0.5, 
       c("CCF","DCF", "ICCF"),
       lty=c(1,1,1),
       lwd=c(2.5,2.5,2.5),col=c("blue","black","red")) 

# ----------------------------------
# insert some NA's

ts.1[10:80,2] <- NA
ts.2[50:100,2] <- NA

plot(ts.1, type="o", bty="n", pch=16, main="chunk of data missing (NA)")
points(ts.2, col="red", type="o")

cat('-- DCF/CCF on data with gaps.', fill=TRUE)
dcf.out <- cross.correlate(ts.1, ts.2, dtau=1, max.lag=50)
dtau <- diff(dcf.out$tau[1:2])
plot(dcf.out$tau-dtau/2, dcf.out$ccf, type="s", lwd=2, bty="n", ylim=c(-1, 1),
     xlab="lag", ylab="CCF", col="black", 
     main="CCF comparison [missing data]\nTrue lag is -10.0")
grid(col="darkgrey")

cat('-- ICCF/CCF on data with gaps.', fill=TRUE)
iccf.out <- cross.correlate(ts.1, ts.2, max.lag=50, dtau=0.5, method="iccf")
lines(iccf.out$tau, iccf.out$ccf, col="red", lwd=2)

ccf.out <- ccf(x, y, plot=FALSE, lag.max=100, na.action = na.pass)
lines(ccf.out$lag, ccf.out$acf, col="blue", lwd=2)

legend(-50, -0.5, 
       c("CCF","DCF", "ICCF"),
       lty=c(1,1,1),
       lwd=c(2.5,2.5,2.5),col=c("blue","black","red")) 

# ----------------------------------
# test ACFs

cat('-- DCF/ACF on data with gaps.', fill=TRUE)
dcf.out <- cross.correlate(ts.1, max.lag=50, method="dcf", dtau=1.0)
dtau <- diff(dcf.out$tau[1:2])
plot(dcf.out$tau-dtau/2, dcf.out$ccf, type="s", lwd=2, bty="n", ylim=c(-1, 1),
     xlab="lag", ylab="CCF", col="black", main="ACF comparison [missing data]")
grid(col="darkgrey")

# now swap ACF and compare
lines(-(dcf.out$tau-dtau/2), dcf.out$ccf, col="pink", lty=3, type="s")

cat('-- ICCF/ACF on data with gaps.', fill=TRUE)
iccf.out <- cross.correlate(ts.1, method="iccf", max.lag=50)
lines(iccf.out$tau, iccf.out$ccf, col="pink", lwd=3)
lines(-iccf.out$tau, iccf.out$ccf, col="black", lwd=2, lty=3)

ccf.out <- ccf(x, x, plot=FALSE, lag.max=100)
lines(ccf.out$lag, ccf.out$acf, col="blue", lwd=2)

legend(-50, -0.5, 
       c("CCF","DCF", "ICCF"),
       lty=c(1,1,1),
       lwd=c(2.5,2.5,2.5),col=c("blue","black","red")) 

# ----------------------------------

# ----------------------------------------
# test CCFs on irregular data

# now remove some data

indx <- sample(1:n, floor(n/2), replace=FALSE)
indx <- sort(indx)
ts.1 <- ts.1[indx,]

indx <- sample(1:n, floor(n/2), replace=FALSE)
indx <- sort(indx)
ts.2 <- ts.2[indx,]

# jitter the times
ts.2$t <- jitter(ts.2$t, amount=0.1)

plot(ts.1, type="o", bty="n", pch=16, main="More data missing\nTime positions jittered")
points(ts.2, col="red", type="o")

dcf.out <- cross.correlate(ts.1, ts.2, max.lag=50, dtau=2.0)
dtau <- diff(dcf.out$tau[1:2])
plot(dcf.out$tau-dtau/2, dcf.out$ccf, type="s", lwd=2, bty="n", ylim=c(-1, 1),
     xlab="lag", ylab="CCF", col="black", 
     main="CCF comparison [missing data]\nTrue lag is -10.0")
grid(col="darkgrey")

iccf.out <- cross.correlate(ts.1, ts.2, max.lag=50,
                            local.est=FALSE, method="iccf")
lines(iccf.out$tau, iccf.out$ccf, col="red", lwd=2)

legend(-50, -0.5, 
       c("DCF", "ICCF"),
       lty=c(1,1),
       lwd=c(2.5,2.5),col=c("black","red")) 

# ----------------------------------
# test ACFs

dcf.out <- cross.correlate(ts.1, max.lag=50, dtau=2.0)
dtau <- diff(dcf.out$tau[1:2])
plot(dcf.out$tau-dtau/2, dcf.out$ccf, type="s", lwd=2, bty="n", ylim=c(-1, 1),
     xlab="lag", ylab="CCF", col="black", main="ACF comparison [missing data]")
grid(col="darkgrey")

# now swap ACF and compare
lines(-(dcf.out$tau-dtau/2), dcf.out$ccf, col="pink", lty=3, type="s")

iccf.out <- cross.correlate(ts.1, ts.1, method="iccf", max.lag=50)
lines(iccf.out$tau, iccf.out$ccf, col="pink", lwd=3)
lines(-iccf.out$tau, iccf.out$ccf, col="black", lwd=2, lty=3)

legend(-50, -0.5, 
       c("DCF", "ICCF"),
       lty=c(1,1),
       lwd=c(2.5,2.5),col=c("black","red")) 

# ----------------------------------
