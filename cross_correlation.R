# -----------------------------------------------------------
cross.correlate <- function(ts.1, ts.2,
                            method = "dcf",
                            max.lag = NULL, 
                            lag.bins = NULL, 
                            min.pts = 5,
                            dtau = NULL,
                            local.est = FALSE,
                            zero.clip = NULL,
                            use.errors = FALSE,
                            one.way = FALSE,
                            cov = FALSE,
                            prob = 0.1, 
                            nsim = 0,
                            peak.frac = 0.8,
                            chatter = 0,
                            plot = FALSE, ...) {
  
  # check arguments
  if (missing(ts.1)) stop('Missing ts.1 data frame.')
  
  # if only one time series then duplicate (compute ACF)
  acf.flag <- FALSE
  if (missing(ts.2)) {
    acf.flag <- TRUE
    ts.2 <- ts.1
  }

  # check the contents of the input time series
  if (!exists("t", where = ts.1)) stop('Missing column t in ts.1.')
  if (!exists("y", where = ts.1)) stop('Missing column y in ts.1.')
  if (!exists("t", where = ts.2)) stop('Missing column t in ts.2.')
  if (!exists("y", where = ts.2)) stop('Missing column y in ts.2.')
  
  # strip out bad data (any NA, NaN or Inf values)
  goodmask.1 <- is.finite(ts.1$t) & is.finite(ts.1$y)
  ts.1 <- ts.1[goodmask.1, ]
  t.1 <- ts.1$t
  n.1 <- length(t.1)
  goodmask.2 <- is.finite(ts.2$t) & is.finite(ts.2$y)
  ts.2 <- ts.2[goodmask.2, ]
  t.2 <- ts.2$t
  n.2 <- length(t.2)
  
  # warning if there's too little data to bother proceeding
  if (n.1 <= 5) stop('ts.1 is too short.')
  if (n.2 <= 5) stop('ts.2 is too short.')

  # if max.tau not set, set it to be 1/4 the max-min time range
  if (is.null(max.lag)) 
    max.lag <- (max(c(t.1, t.2)) - min(c(t.1, t.2))) * 0.25
  
  # if dtau is not defined, set to default
  if (is.null(dtau))
    dtau <- min( diff(t.1) )
  
  # total number of lag bins; make sure its odd!
  n.tau <- ceiling(2 * max.lag / dtau)  
  if (n.tau %% 2 == 0) 
    n.tau <- n.tau + 1
  lag.bins <- round( (n.tau - 1) / 2 )
  if (lag.bins < 4) stop('You have too few lag bins.')
  
  # now adjust max.lag so that -max.lag to +max.lag spans odd number of 
  # dtau bins
  max.lag <- (1/2) * dtau * (n.tau-1)
  
  # define the vector of lag bins. tau is the centre of each bin.
  # should extend from -max.tau to +max.tau, centred on zero.
  tau <- seq(-max.lag, max.lag, by = dtau)
  
  # optional feedback for the user  
  if (chatter > 0) {
    cat('-- lag.bins:', lag.bins, fill = TRUE)
    cat('-- max.lag:', max.lag, fill = TRUE)
    cat('-- n.tau:', n.tau, fill = TRUE)
    cat('-- length(tau)', length(tau), fill = TRUE)
    cat('-- dtau:', dtau, fill = TRUE)
    cat('-- length ts.1:', n.1, ' dt:', diff(t.1[1:2]), fill = TRUE)
    cat('-- length ts.2:', n.2, ' dt:', diff(t.2[1:2]), fill = TRUE)
  }
  
  # compute CCF for the input data
  ccf.out <- NA
  if (method == "dcf") {
    ccf.out <- dcf(ts.1, ts.2, tau,
                   local.est = local.est, min.pts = min.pts, 
                   zero.clip = zero.clip, chatter = chatter, cov = cov)
  } else {
    ccf.out <- iccf(ts.1, ts.2, tau,
                    local.est = local.est, chatter = chatter, 
                    cov = cov, one.way = one.way, zero.clip = zero.clip)
  }
  
  # extract the lag settings
  max.lag <- max(ccf.out$tau)
  nlag <- length(ccf.out$tau)
  lag.bins <- (nlag-1)/2
  dtau <- diff(ccf.out$tau[1:2])
  
  if (chatter > 0) {
    cat('-- max.lag:   ', max.lag, fill = TRUE)
    cat('-- nlag:      ', nlag, fill = TRUE)
    cat('-- lag.bins:  ', lag.bins, fill = TRUE)
    cat('-- dtau:      ', dtau, fill = TRUE)
  }
  
  # plot the CCF
  if (plot == TRUE) {
    plot(0, 0, type = "n", lwd = 2, bty = "n", ylim = c(-1, 1), xlim = c(-1, 1)*max.lag,
         xlab = "lag", ylab = "CCF", ...)
    grid(col = "lightgrey")
    dtau <- diff(ccf.out$tau[1:2])
    if (method == "dcf") {
      lines(ccf.out$tau-dtau/2, ccf.out$ccf, type="s", lwd = 2, col = "blue")
    } else {
      lines(ccf.out$tau, ccf.out$ccf, col = "blue", lwd = 3)
    }
  }
  
  # set up blanks if no simulations are run
  lower <- NA
  upper <- NA
  cent.dist <- NA
  peak.dist <- NA
  
  # (optional) run simulations to compute errors
  if (nsim > 0) {
  
    # check we have enough simulations
    if (nsim < 2/max(prob)) stop('Not enough simulations. Make nsim or prob larger.')
        
    # run some simulations
    sims <- ccf.errors(ts.1, ts.2, tau, nsim = nsim,
                     method = method, peak.frac = peak.frac, min.pts = min.pts,
                     local.est = local.est, zero.clip = zero.clip, prob = prob,
                     cov = cov, chatter = chatter, acf.flag = acf.flag, 
                     one.way = one.way)
    
    lower <- sims$lags$lower
    upper <- sims$lags$upper
    cent.dist <- sims$dists$cent.lag
    peak.disk <- sims$dists$peak.lag
  
  } else {
    
    if (method == "dcf") {
      acf.1 <- dcf(ts.1,  ts.1, tau,
                   local.est = local.est, 
                   min.pts = min.pts, 
                   chatter = chatter, cov = cov)
      acf.2 <- dcf(ts.2, ts.2,  
                   local.est = local.est, 
                   min.pts = min.pts, dtau = dtau,
                   max.lag = max.lag, lag.bins = lag.bins, 
                   chatter = chatter, cov = cov)
    } else {
      acf.1 <- iccf(ts.1, ts.1, tau, 
                    local.est = local.est, 
                    chatter = chatter, 
                    cov = cov, 
                    one.way=one.way)
      acf.2 <- iccf(ts.2, ts.2, tau, 
                    local.est = local.est, 
                    chatter = chatter, 
                    cov = cov, 
                    one.way = one.way)
    }
    
    sigma <- sqrt( (1/ccf.out$n) * sum(acf.1$ccf*acf.2$ccf) )
    lower <- ccf.out$ccf - sigma
    upper <- ccf.out$ccf + sigma
  }
  
  # plot confidence bands
  if (plot == TRUE) {
    pnk <- rgb(255, 192, 203, 100, maxColorValue = 255)
    indx <- is.finite(ccf.out$ccf)
    polygon(c(ccf.out$tau[indx], rev(ccf.out$tau[indx])), c(lower[indx], rev(upper[indx])), 
            col=pnk, border = NA)
    if (method == "dcf") {
      lines(ccf.out$tau-dtau/2, ccf.out$ccf, type = "s", lwd = 2, col = "blue")
    } else {
      lines(ccf.out$tau, ccf.out$ccf, col = "blue", lwd = 3)
    }
  }
  
  # return output  
  result <- list(tau = ccf.out$tau, 
                 ccf = ccf.out$ccf, 
                 lower = lower,
                 upper = upper, 
                 peak.dist = peak.dist,
                 cent.dist = cent.dist, 
                 method = method)
  return(result)
}

# -----------------------------------------------------------
# matrix.tau
# Inputs: 
#   t.1   - times for vector x.1
#   t.2   - times for vector x.2
#
# Value:
#  tau.ij - matrix of time lags 
#
# Description:
#  Define the matrix of time lags
#   tau[i,j] = t.1[i] - t.2[j]
#
# Note that in the special case that t.1=t.2 we have
# a square symmetric matrix: tau[i,j] = tau[j,i]. 
# In the even more special case that
# the two series are identically and evenly sampled
# (t.1[i] = t.2[i] = i * dt + t0) then we have a 
# circulant matrix; the jth column tau[,j] is the (j-1)th 
# cyclic permutation of the first column. This matrix is
# symmetric, Toeplitz and circulant. 
#
# History:
#  21/03/16 - First working version
#
# (c) Simon Vaughan, University of Leicester
# -----------------------------------------------------------

matrix.tau <- function(t.1, t.2) {

  # check arguments
  if (missing(t.1)) stop('Missing t.1 vector.')
  if (missing(t.2)) t.2 <- t.1 
  
  # compute t.0, an arbitrary start time
  # this improves accuracy if the time offset is large
  # compared to the range of times.
  t.0 <- min(t.1, t.2)
  t.1 <- t.1 - t.0
  t.2 <- t.2 - t.0
  
  # compute the time lag matrix
  tau <- outer(t.1, t.2, "-")
  
  # return to calling function
  return(tau)  
}

# -----------------------------------------------------------
# dcf
# Inputs: 
#   ts.1      - data frame containing times (t)
#                and values (x) for data series 1.
#                Contains optional errors (dx).
#   ts.2      - data frame for data series 2
#   tau       - list of lag bins
#   min.pts   - set to NA any lag bins with fewer points (default: 10)
#   local.est - use 'local' (not 'global') means and variances?
#   zero.clip - remove pairs of points with exactly zero lag? 
#   use.errors - TRUE/FALSE: if TRUE then subtract mean square error from variances
#   cov       - compute covariances rather than correlations (default: FALSE)
#   chatter   - (integer) level of information reported while running
#
# Value:
#   result    - a data frame containing columns...
#      tau    - the centre of the lag bins (vector)
#      ccf    - the correlation coefficent in each lag bin
#
# Description:
#  Compute the Discrete Correlation Function based on the method 
# outlined in Edelson & Korlik (1998, ApJ, v333, pp646-659). 
#
# This is a way to estimate the CCF for a pair of time series
# (t.1, x.1) and (t.2, x.2) when the time sampling is uneven
# and non-synchronous. 
#
# Input are two time series (data frames with columns: t, x, dx [optional])
# Output is the correlation coefficient r(tau) in different lag bins.
# A peak in t at positive lag indicates that ts.2 leads ts.1; a negative lag 
# indicates that ts.1 leads ts.2.
#
# We first subtract the mean values from x.1 and x.2. Then,
# within the ith lag bin (tau[i] - dtau/2, tau[i] + dtau/2) we collect 
# all pairs (x.1, x.2) of data for which t.1 - t.2 falls within the lag bin.
# Using these pairs of data we compute the sum of their product and
# normalise it:
#
#    cov[i] = (1/n) * sum_k x.1[k] * x.2[k] 
#
# Here, k is the set of values of points within the lag bin and has
# n pairs. This gives a covariance. If cov=TRUE we keep these values.
# Otherwise (default: cov=FALSE) we normalise by the product of the 
# standard deviations of x.1 and x.2: sqrt(var(x.1)*var(x.2)).
#
# The number of output lags (tau) is 2*lag.bins+1, and the lags are
# centred on zero. So they run from -max.lag to +max.lag.
# Any lag bins containing fewer than min.pts pairs of points will have
# dcf = NA. IF lag.bins is not supplied it is set to ~1/4 the number of
# points in the short time series.
#
# If max.lag is not supplied it will be set to 1/4 the difference between 
# the very first and last times (from either time series).
#
# If local.est == FALSE (default) then the correlation coefficient is computed
# using the 'global' mean and variance of each time series. 
#
# If local.est == TRUE then the correlation coefficient is computed
# using the 'local' mean and variance. Within each lag bin the 
# mean to be subtrated and the varaince to be divided are computed using
# only data points contributing to that lag bin, i.e. using only x.1 and x.2 for 
# which the corresponding t.1 - t.2 is within the range 
# (tau[i] - dtau/2, tau[i] + dtau/2).
#
# If 'errors' are suppled for either or both time series (dx.1, dx.2) 
# and you have set use.errors == TRUE then the variance used 
# in the denominator terms of correlation coefficient calculation will be
# the 'excess variance', i.e. the total sample variance minus the 
# mean square error. If dx.1 and/or dx.2 are a single number, this is 
# assumed to be the same error for each data point. If use.errors == FALSE
# (default) then the usual sample variance will be used.
#
# If zero.clip=TRUE then pairs of points with zero lag between them
# are not included in the correlation estimates. This helps remove
# zero-lag errors, when measurements in two bands are affected by systematic
# errors occuring at the same time.
#
# History:
#  21/03/16 - v1.0 - First working version
#  05/04/16 - v1.1 - added na.rm option to strip out non-finite values
#  09/04/16 - v1.2 - added use.errors option; minor fixes
#  11/04/16 - v1.3 - added dtau and chatter input; bug fixes;
#                     removed na.rm (now automatic)
#  24/07/16 - v1.4 - moved most checks and tau calculation to the 
#                     main cross.correlation function.
#
# (c) Simon Vaughan, University of Leicester
# -----------------------------------------------------------

dcf <- function(ts.1, ts.2, 
                tau= NULL,
                min.pts = 5,
                local.est = FALSE,
                zero.clip = NULL,
                use.errors = FALSE,
                cov = FALSE,
                chatter = 0) {

  # check arguments
  if (missing(ts.1)) stop('Missing ts.1 data frame.')
  if (missing(ts.2)) stop('Missing ts.2 data frame')
  if (is.null(tau)) stop('Missing tau in.')
  
  # extract columns from ts.1
  # create "error" column if not present
  t.1 <- ts.1$t
  x.1 <- ts.1$y
  n.1 <- length(t.1)
  if (exists("dy", where = ts.1)) {
    dx.1 <- ts.1$dy
  } else {
    dx.1 <- rep(0, n.1)
  }
  
  # extract columns from ts.2
  # create "error" column if not present
  t.2 <- ts.2$t
  x.2 <- ts.2$y
  n.2 <- length(t.2)
  if (exists("dy", where = ts.2)) {
    dx.2 <- ts.2$dy
  } else {
    dx.2 <- rep(0, n.2)
  }
  
  if (length(dx.1) == 1) dx.1 <- rep(dx.1, n.1)
  if (length(dx.2) == 1) dx.2 <- rep(dx.2, n.2)
  
  # remove any (possible) large offset on the time stamps
  # this may improve numerical accuracy on lags if the t.0
  # is very large compared to the sampling interval.
  t.0 <- min(c(t.1, t.2))
  t.1 <- t.1 - t.0
  t.2 <- t.2 - t.0
  
  # subtract the global mean from each time series
  x.1 <- x.1 - mean(x.1)
  x.2 <- x.2 - mean(x.2)
  
  # compute variance of each series
  var.1 <- var(x.1, na.rm = TRUE)
  var.2 <- var(x.2, na.rm = TRUE)
  
  # set mean square error to zero unless use.errors == TRUE
  # else compute mean square error and remove from variance
  if (use.errors == FALSE) {
    sd.1 <- sqrt( var.1 )
    sd.2 <- sqrt( var.2 )
  } else {
    mse.1 <- mean( dx.1^2 )
    mse.2 <- mean( dx.2^2 )
    sd.1 <- sqrt( var.1 - mse.1 )
    sd.2 <- sqrt( var.2 - mse.2 )
  }
  
  # define lag settings
  n.tau <- length(tau)
  lag.bins <- as.integer(n.tau - 1)/2
  dtau <- diff(tau[1:2])
  
  # compute a n.1 x n.2 matrix of the lags for every pair of data points
  tau.12 <- matrix.tau(t.1, t.2)
  
  # convert matrix of tau values to lag bins (-lag.bins:lag.bins)
  # by rounding the lag to the nearest integer of dtau.
  tau.12 <- round( tau.12 / dtau ) 
  
  # The matrix of lags is forced to by of integer type to save 
  # memory, and to speed up the slow subset selection step later.
  storage.mode(tau.12) <- "integer"
  
  # prepare the vector to store the correlation coefficient values
  dcf <- array(NA, dim = n.tau)
  n.i <- array(1, dim = n.tau)
  
  # loop through each lag bin
  for (i in 1:n.tau) {
    
    # select only pairs of points whose time differences fall
    # within the range of lag bin i. This is forced to be of 
    # integer type to run faster when comparing to the integer type matrix.
    lag.i <- as.integer(round(tau[i] / dtau))
    
    # This is the *slow* part. Suggestions welcome....
    indx <- which(tau.12 == lag.i, arr.ind = TRUE)

    # extract the relevant data points
    x1.i <- x.1[indx[, 1]]
    x2.i <- x.2[indx[, 2]]
    dx1.i <- dx.1[indx[, 1]]
    dx2.i <- dx.2[indx[, 2]]
    n.i[i] <- length(x1.i)

    # count the number of pairs are zero lag
    if ( i == lag.bins+1 ) n.0 <- n.i[i]
    
    # if using zero clipping, then also ignore all pairs with t_i - t_j = 0
    if ( !is.null(zero.clip) & i == lag.bins+1 ) {
      #indx <- which(tau.12 == lag.i & tau.12 != 0, arr.ind = TRUE)
      t1.i <- t.1[indx[, 1]]
      t2.i <- t.2[indx[, 2]]
      mask <- ( abs(t1.i - t2.i) > dtau/2*zero.clip)
      x1.i <- x1.i[mask]
      x2.i <- x2.i[mask]
      dx1.i <- dx1.i[mask]
      dx2.i <- dx2.i[mask]
      n.i[i] <- length(x1.i)
    }
    
    # skip any lag bin with fewer than min.pts pairs of points
    if (n.i[i] < min.pts) next
    
    # subtract means (local or global) from the data within lag bin i
    if (local.est == TRUE & cov == FALSE) {
      mean.1 <- mean(x1.i)
      mean.2 <- mean(x2.i)
      x1.i <- x1.i - mean.1
      x2.i <- x2.i - mean.2
      var.1 <- var(x1.i)
      var.2 <- var(x2.i)
      if (use.errors == FALSE) {
        sd.1 <- sqrt( var.1  )
        sd.2 <- sqrt( var.2  )
      } else {
        mse.1 <- mean( dx1.i^2 )
        mse.2 <- mean( dx2.i^2 )
        sd.1 <- sqrt( var.1 - mse.1 )
        sd.2 <- sqrt( var.2 - mse.2 )
      }
    }

    # if cov == TRUE then compute covariance, not correlation
    # I.e. do not normalise by the sqrt(variances)
    if (cov == TRUE) {
      sd.1 <- 1.0
      sd.2 <- 1.0
    }
    
    # now compute the product moment
    dcf[i] <- sum( x1.i * x2.i ) / (sd.1 * sd.2) 
  }

  dcf <- dcf / n.i

  # health warning
  if (chatter > 0) {
    r.max <- max(abs(dcf), na.rm = TRUE)
    if (r.max > 1 & cov == FALSE)
      warning(paste('CCF outside of (-1,+1): ', r.max))
  }

  # return final output
  return(data.frame(tau = tau, ccf = dcf, n = n.i))
  
}

# --------------------------------------------------
# NAME: 
#     iccf
#
# PURPOSE:
#     Estimate a cross correlation between two time series
#
# AUTHOR:
#     Simon Vaughan
#
# CALLING SEQUENCE:
#     r <- iccf(ts.1, ts.2, tau)
#
# INPUTS:
#   ts.1      - (data frame) times series 1, contains t, y.
#   ts.2      - (data frame) times series 2, contains t, y.
#   tau       - list of lags
#   cov       - compute covariances rather than correlations (default: FALSE)
#   chatter   - (integer) level of information reported while running
#
# Value:
#   result    - a data frame containing columns...
#      tau    - the centre of the lag bins (vector)
#      ccf    - the correlation coefficent in each lag bin
#
# DETAILS:
#
# Given two time series x.1 and x.2 sampled at time t.1 and t.2
# we estimate a cross correlation function (CCF) by interpolating
# (t.2, y.2). For a lag tau we estimate y.2 at each time t.1+tau
# by interpolating between the two nearest points of y.2. We then
# pair the values y.1 with the corresponding lagged values of y.2
# and compute the linear correlation coefficient, r. We repeat this
# for a range of lags and plot r vs. tau.
#
# The interpolation is handled by the approx() function
# for linear interpolation.
#
# History:
#  21/03/16 - v1.0 - First working version
#  05/04/16 - v1.1 - added na.rm option to strip out non-finite values
#  09/04/16 - v1.2 - minor fixes
#  24/07/16 - v1.3 - moved most checks and tau calculation to the 
#                     main cross.correlation function.
#
# (c) Simon Vaughan, University of Leicester
# -----------------------------------------------------------

iccf <- function(ts.1, ts.2, 
                 tau = NULL,
                 local.est = FALSE,
                 one.way = FALSE,
                 zero.clip = NULL,
                 cov = FALSE,
                 chatter = 0) {
  
  # check arguments
  if (missing(ts.1)) stop('Missing TS.1 argument ICCF')
  if (missing(ts.2)) stop('Missing TS.2 argument ICCF')
  if (is.null(tau)) stop('Missing tau in.')
  
  t.1 <- ts.1$t
  x.1 <- ts.1$y
  n.1 <- length(t.1)

  t.2 <- ts.2$t
  x.2 <- ts.2$y
  n.2 <- length(t.2)

  # remove means
  x.1 <- x.1 - mean(x.1)
  x.2 <- x.2 - mean(x.2)

  # define lag settings
  n.tau <- length(tau)
  lag.bins <- as.integer(n.tau - 1)/2
  dtau <- diff(tau[1:2])

  # main calculation
  iccf.ij <- iccf.core(t.1, x.1, t.2, x.2, tau, local.est = local.est, cov = cov)
  r.ij <- iccf.ij$r
  if (one.way == FALSE) {
    iccf.ji <- iccf.core(t.2, x.2, t.1, x.1, -tau, local.est = local.est, cov = cov)
    r.ji <- iccf.ji$r
    r <- 0.5 * (r.ij + r.ji)
  } else {
    r <- r.ij
  }
  
  # health warning
  if (chatter > 0) {
    r.max <- max(abs(r), na.rm = TRUE)
    if (r.max > 1 & cov == FALSE)
      warning(paste('CCF outside of (-1, +1): ', r.max))
  }
  
  # remove the zero-lag bin if requested
  if (!is.null(zero.clip))
    r[lag.bins+1] <- NA
  
  # return output to user  
  return(data.frame(tau = tau, ccf = r, n = iccf.ij$n))
  
}

# -----------------------------------------------------------
# The main loop for the ICCF
# In this part we take time series 1, x.1 at t.1, pair them with values
# from time series 2, x.2 at t.1-tau[i] produce by linearly interpolating
# between the nearest values of x.2.
# At a given tau[i] we su, the product of the paired x.1 and x.2 values 
#  r[i] = (1/n) * sum(x.1 * x.2) / (sd.1 * sd.2)
# In the simplest case n, sd.1 and sd.2 are constant and are the
# number of pairs at lag=0 and the total sqrt(var) of each time series.
# If local.est is TRUE then n, sd.1 and sd.2 are evaluated "locally"
# i.e. they are vary for each lag tau[i]. In this case they are the
# number of good pairs at lag tau[i], and the sqrt(vars) of just the 
# x.1 and x.2 data points involved.
# We assume x.1 and x.2 have zero sample mean.

iccf.core <- function(t.1, x.1, 
                      t.2, x.2, 
                      tau, 
                      local.est = FALSE,
                      cov = FALSE) {
  
  n.tau <- length(tau)
  lag.bins <- as.integer( (n.tau-1)/2 )
  
  r.ij <- array(0, dim = n.tau)     # correlation of x1(t_i) vs. x2(t_2)
  n.ij <- array(1, dim = n.tau)     # no. data pairs at lag t_i - t_j

  sd.1 <- sd(x.1)                 # sqrt(var) of time series 1
  sd.2 <- sd(x.2)                 # sqrt(var) of time series 2
  
  for (i in 1:n.tau) {
    
    # at the ith lag, tau[i], estimate the values of the second time series
    # x.2 at times t.1 - tau[i] by interpolation between nearby points.
    # Note, values interpolated outside of the range of x.2 are set to NA.
    x.2interp <- approx(t.2, x.2, t.1-tau[i])$y

    # select only complete pairs of data (remove NA's)
    indx.ij <- is.finite(x.2interp) & is.finite(x.1)
    xi.1 <- x.1[indx.ij]
    xi.2 <- x.2interp[indx.ij]

    # compute the local mean and variance if needed
    if (local.est == TRUE & cov == FALSE) {
      xi.1 <- xi.1 - mean(xi.1)
      xi.2 <- xi.2 - mean(xi.2)
      sd.1 <- sd(xi.1)
      sd.2 <- sd(xi.2)
    }
 
    # if cov == TRUE then compute covariance, not correlation
    # I.e. do not normalise by the sqrt(variances)
    if (cov == TRUE) {
      sd.1 <- 1.0
      sd.2 <- 1.0
    }
    
    # compute sum of the products of good pairs
    n.ij[i] <- length(xi.1)
    r.ij[i] <- sum(xi.1 * xi.2) / sd.1 / sd.2

  }
  
  # if local.est == FALSE then normalise by pairs at lag = 0.
  # if local.est == TRUE use actual no pairs in lag bin i.
  if (local.est == TRUE) {
    r.ij <- r.ij / n.ij
  } else {
    r.ij <- r.ij / n.ij[lag.bins+1]
  }
  
  # return vector of coefficients 
  return(list(r=r.ij, n=n.ij))

}

# -----------------------------------------------------------
# fr.rss
# Inputs: 
#   ts.1      - data frame containing times (1st column)
#                and values (2nd column) for data series.
#                Contains optional errors (3rd column).
#
# Value:
#   result   - a data frame containing columns
#      t     - time bins for randomised data
#      x     - values for randomised data
#     dx     - errors for randomised data
#
# Description:
#  Performs "flux randomisation" and "random sample selection"
# of an input time series, following 
# Peterson et al. (2004, ApJ, 613:682-699).
#
# Given an input data series (t, x, dx) of length N we sample
# N points with replacement. Duplicated points are ignored, so 
# the ouptut is usually shorter than the input. So far this is
# a basic bootstrap procedure.
#
# If error bars are provided: when a point is selected m times, 
# we decrease the error by 1/sqrt(m). See Appendix A of Peterson et al. 
# And after resampling in time, we then add a random Gaussian deviate
# to each remaining data point, with std.dev equal to its error bar.
# In this way both the times and values are randomised.
#
# If errors bars are nor provided, this is a simple bootstrap.
#
# The output is another data frame of (t, y, dy)
#
# History:
#  21/03/16 - First working version
#
# (c) Simon Vaughan, University of Leicester
# -----------------------------------------------------------

fr.rss <- function(dat) {
  
  # check arguments
  if (missing(dat)) stop('Missing DAT argument')
  
  # extract data
  times <- dat$t
  x <- dat$y
  n <- length(times)
  
  # if no errors are provided, we will do a simple bootstrat
  if (exists("dy", where = dat)) {
    bootstrap <- FALSE
    dx <- dat$dy
  } else {
    bootstrap <- TRUE
    dx <- 0
  }
  if (length(dx) < n) dx <- rep(dx[1], n)
    
  # ignore any data with NA value
  mask <- is.finite(x) & is.finite(dx)
  times <- times[mask]   # time
  x <- x[mask]           # value
  dx <- dx[mask]         # error
  n <- length(x)

  # randomly sample the data
  indx <- sample(1:n, n, replace=TRUE)

  # identify which points are sampled more than once
  duplicates <- duplicated(indx)
  indx.clean <- indx[!duplicates]

  # where data points are selected n>0 times, scale the error by 1/sqrt(n)
  n.repeats <- hist(indx, plot=FALSE, breaks=0:n+0.5)$counts
  dx.original <- dx
  dx <- dx / sqrt( pmax.int(1, n.repeats) )
    
  # new data are free from duplicates, and have errors decreased where
  # points are selected multiple times.
  
  times.new <-  times[indx.clean]
  x.new <-      x[indx.clean]
  dx.new <-     dx[indx.clean]
  n <- length(x.new)
  
  # now randomise the fluxes at each data point
  if (bootstrap == FALSE) {
    x.new <- rnorm(n, mean=x.new, sd=dx.new)
  }      
  
  # sort into time order
  indx <- order(times.new)

  # return the ouput data
  return(data.frame(t=times.new[indx], y=x.new[indx], dy=dx.new[indx]))
  
}

# -----------------------------------------------------------
# ccf.errors
# Inputs: 
#   ts.1      - data frame containing times (t)
#                and values (y) for data series 1.
#                Contains optional errors (3rd column).
#   ts.2      - data frame for data series 2
#   tau       - list of lags
#   min.pts   - set to NA any lag bins with fewer points (default: 10)
#   local.est - use 'local' (not 'global') means and variances?
#   prob      - probability level to use for confidence intervals
#   nsim      - number of simulations
#   peak.frac - what fraction below peak to include in centroid measurements?
#   zero.clip - remove pairs of points with exactly zero lag? 
#   method    - use DCF or ICCF method (default: dcf)
#   use.errors - TRUE/FALSE passed to dcf() 
#   local.est  - TRUE/FALSE passed to dcf() or iccf()
#   acf.flag  - TRUE if computing ACF and ts.2 = ts.1
#
# Value:
#   result    - a data frame containing columns...
#      tau    - the centre of the lag bins (vector)
#      dcf    - the correlation coefficent in each lag bin
#
# Description:
#  Compute the Discrete Correlation Function based on method 
# outlined in Edelson & Korlik (1998, ApJ). 
# This is a way to estimate the CCF for a pair of time series
# (t.1, x.1) and (t.2, x.2) when the time sampling is uneven
# and non-synchronous. See dcf(...) function.
#
# Computes errors on the DCF values using "flux randomisation" 
# and "random subset sampling" FR/RSS using the fr.rss(...) function.
#
# For each randomised pair of light curves we compute the DCF. 
# We record the DCF, the lag at the peak, and the centroid lag
# (including only points higher than peak.frac * peak).
# Using nsim simulations we compute the (1-p)*100% confidence
# intervals on the DCF values, and the distribution of peaks
# and centroids.
#
# The output is a list containing two components
#    lags     - a data frame with four columns
#    tau      - time lags
#    dcf      - the DCF values for the input data
#    lower    - the lower limit of the confidence interval
#    upper    - the upper limit of the confidence interval
#    dists    - a data frame with two columns
#    peak.lag - the peak values from nsim simulations
#    cent.lag - the centroid values from nsim simulations
#
# History:
#  21/03/16 - v1.0 - First working version
#  05/04/16 - v1.1 - added na.rm option to strip out non-finite values
#  09/04/16 - v1.2 - added use.errors option; minor fixes
#  23/07/16 - v1.3 - minor change to handling of centroid calculation.
#                     if the CCF is entirely <0 then return NA for
#                     centroid.
#
# (c) Simon Vaughan, University of Leicester
# -----------------------------------------------------------

ccf.errors <- function(ts.1, ts.2, 
                       tau = NULL,
                       min.pts=5,
                       local.est=FALSE,
                       cov=FALSE,
                       prob=0.1, 
                       nsim=250,
                       peak.frac=0.8,
                       zero.clip=NULL,
                       one.way=FALSE,
                       method="dcf",
                       use.errors=FALSE,
                       acf.flag=FALSE,
                       chatter=0) {

  # check arguments
  if (missing(ts.1)) stop('Missing ts.1 data frame.')
  if (is.null(tau)) stop('Missing tau in.')
  
  # if only one time series then duplicate (compute ACF)
  acf.flag <- FALSE
  if (missing(ts.2)) {
    acf.flag <- TRUE
    ts.2 <- ts.1
  }

  if (peak.frac > 1 | peak.frac < 0) 
    stop('peak.frac should be in the range 0-1')
  
  if (prob > 1 | prob < 0) 
    stop('prob should be in the range 0-1')
  
  nlag <- length(tau)

  # set up an array for the simulated DCFs
  ccf.sim <- array(NA, dim=c(nsim, nlag))  
  peak.lag <- array(NA, dim=nsim)
  cent.lag <- array(NA, dim=nsim)

  # loop over simulations
  for (i in 1:nsim) {
    
    # generate randomised data
    ts1.sim <- fr.rss(ts.1)
    if (acf.flag == FALSE) {
      ts2.sim <- fr.rss(ts.2)
    } else {
      ts2.sim <- ts1.sim
    }
    
    # compute CCF of randomised data
    if (method == "dcf") {
      result.sim <- dcf(ts1.sim, ts2.sim, tau,
                        local.est = local.est, 
                        min.pts = min.pts, 
                        cov = cov,
                        zero.clip = zero.clip, 
                        use.errors = use.errors)
    } else {
      result.sim <- iccf(ts1.sim, ts2.sim, tau,
                         zero.clip = zero.clip, 
                         local.est = local.est, 
                         cov = cov, 
                         one.way = one.way)
    }
    
    ccf.sim[i,] <- result.sim$ccf    
    
    # find and store the peak of the CCF
    peak <- which.max(result.sim$ccf)
    peak.lag[i] <- result.sim$tau[peak]
    
    # find and store the centroid of the CCF
    # if the CCF peak is <0 then return NA 
    if (max(result.sim$ccf, na.rm = TRUE) > 0) {
      mask <- which( result.sim$ccf >= peak.frac * max(result.sim$ccf, na.rm = TRUE) )
      cent.lag[i] <- sum(result.sim$ccf[mask]*result.sim$tau[mask], na.rm = TRUE)  /
                      sum(result.sim$ccf[mask], na.rm = TRUE)
    }
    
    cat("\r Processed", i, "of", nsim, "simulations")
  }

  cat(fill = TRUE)

  # now compute the prob/2 and 1-prob/2 quantiles at each lag
  ccf.lims <- array(NA, dim = c(nlag, 2))
  probs <- c(prob/2, 1-prob/2)
  for (i in 1:nlag) {
    ccf.lims[i,] <- quantile(ccf.sim[,i], probs = probs, na.rm = TRUE)
  }

  # package the DCF: lag, values, lower and upper limits  
  lags <- data.frame(tau = result.sim$tau, lower = ccf.lims[,1], upper = ccf.lims[,2])
  
  # package the peak and centroid data
  dists <- data.frame(peak.lag = peak.lag, cent.lag = cent.lag)
  
  # interval of centroids
  if (chatter >= 0) 
  cat('-- ', signif(100*(1-prob), 3), '% lag interval ', 
      quantile(cent.lag, probs[1], na.rm = TRUE), ' - ',
      quantile(cent.lag, probs[2], na.rm = TRUE), fill = TRUE, sep = "")
  
  
  return(list(lags = lags, dists = dists))
}

# -----------------------------------------------------------
