# define function for ACV
acv <- function(tau, A, tau.0) {
  a <- A*exp(-abs(tau)/tau.0)
  return(a)
}

# define function for power spectrum
psd <- function(f, B, f.0) {
  s <- B/(1+(f/f.0)^2)
  return(s)
}

# define ACV parameters
A <- 1.0
tau.0 <- 3.0

# convert to corresponding B, f.0
f.0 <- 1/(2*pi*tau.0)
B <- 2*A*tau.0

# define grid of tau values
tau <- seq(0, 50, by = 0.5)
n <- length(tau)

# compute the true ACV from 0 to +100
acv.true <- acv(tau, A, tau.0)

# plot the true curve
par(oma=c(0,0,0,0), mar=c(5,4,1,1))
layout(c(1,2), heights=c(0.7,0.3))
max.t <- 50
plot(tau, acv.true, type="l", xlim=c(0,max.t))

# compute the PSD for positive frequencies
# so f = f_0, f_1, ..., f_n-1, f_n
f <- seq(0, 20, by = 0.001)
n <- length(f)

# now add on the negatative frequencies, 
# by 'mirroring' around the top element.
# and then removing the last element.
# So f = f_0, f_1, ..., f_n-1, f_n, -f_n-1, ..., -f_2, -f_1.
# Note that since S(-f) = S(f) we can use f_j as -f_j.
f <- c(f, rev(f))
f <- f[-length(f)]

# resolution of frequency grid
df <- diff(f[1:2])

# compute S(f) for +ve and -ve frequencies
S <- psd(f, B, f.0)

# compute the FFT
result <- fft(S) * df
acv.approx <- Re(result)

# extract only values corresponding to +ve tau's
acv.approx <- acv.approx[1:n]

# compute the corresponding tau's
tau.approx <- (0:(n-1)) / length(f) / df # (2*max(f))

# plot approximation
lines(tau.approx, acv.approx, col="red")

# now compute the true ACV at these tau's
acv.true <- acv(tau.approx, A, tau.0)

# differences
d <- (acv.true - acv.approx) 

# plot the result
plot(tau.approx, d, col="red", xlim=c(0, max.t), 
     ylim=c(-1,1)*max(abs(d)), type="l", xlab="tau", 
     ylab="difference")
