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
