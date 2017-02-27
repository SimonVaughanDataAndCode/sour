# --------------------------------------------------
# iccf
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

#' Estimate a cross correlation between two time series
#' 
#' \code{iccf} returns the Interpolated Cross-Correlation Function estimates.
#'
#' @param tau (vector) list of lags at which to compute the CCF.
#' @inheritParams cross.correlate
#'
#' @return
#' A data frame containing columns:
#'  \item{tau}{lags (in time units)}
#'  \item{ccf}{correlations coefficent in each lag bin}
#'
#' @section Notes:
#' In what follows we refer to the \code{t, y} values of time series 1 
#' (\code{ts.1}) as \code{t.1, y.1}, and similarly for time series 2.
#' 
#' Given two time series \code{y.1} and \code{y.2}, sampled at times \code{t.1}
#' and \code{t.2}, we estimate a cross correlation function (CCF) by 
#' interpolating (\code{t.2, y.2}). For a lag \code{tau} we estimate \code{y.2} 
#' at each time \code{t.1+tau} by interpolating between the two nearest points 
#' of \code{y.2}. We then pair the values \code{y.1} with the corresponding 
#' lagged values of \code{y.2} and compute the linear correlation coefficient, 
#' \code{ccf}. The interpolation is handled by the \code{approx} function for
#' linear interpolation.
#' 
#' @seealso \code{\link{cross.correlate}}, \code{\link{dcf}}, \code{\link[stats]{approx}}
#' 
#' @examples 
#' ## Example using NGC 5548 data
#' res <- iccf(cont, hbeta, tau = seq(-100, 100))
#' plot(res$tau, res$ccf, type = "l", col = "blue", lwd = 3, bty = "n")
#' grid()
#' 
#' @export
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
  
  # main calculation (using iccf.core function)
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

#' Compute the one-way Interpolated Cross-Correlation Function (ICCF)

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
    r.ij <- r.ij / (n.ij-1)
  } else {
    r.ij <- r.ij / (n.ij[lag.bins+1] - 1)
  }
  
  # return vector of coefficients 
  return(list(r=r.ij, n=n.ij))
  
}
