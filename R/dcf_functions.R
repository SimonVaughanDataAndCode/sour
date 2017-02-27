# -----------------------------------------------------------
# matrix.tau
#
# History:
#  21/03/16 - v1.0 - First working version
#  27/02/17 - v1.1 - renamed input variables.
#
# (c) Simon Vaughan, University of Leicester
# -----------------------------------------------------------

#' Compute a matrix of differences given two vectors.
#' 
#' \code{matrix.tau} returns a matrix of differences between two vectors.
#' 
#' Given two vectors - \code{x} (length \code{M}) and \code{y} (length \code{N}) - as 
#' input, return the \code{N*M} matrix of differences \code{result[i,j] = x[i] - y[j]}.
#' 
#' @param x vector 1
#' @param y vector 2 (default to vector 1 if not specified)
#'
#' @return 
#' \code{N*M} array of differences, \code{result[i,j] = x[i] - y[j]}
#'
#' @section Notes:
#' Note that in the special case that \code{x=y} we have a square symmetric
#' matrix: \code{result[i,j] = result[j,i]}. In the even more special case that 
#' the two vectors are evenly spaced (\code{x[i] = y[i] = i * delta + const})
#' then we have a circulant matrix; the \code{j}th column \code{result[,j]} is
#' the \code{(j-1)}th cyclic permutation of the first column. This matrix is 
#' symmetric, Toeplitz and circulant.
#' 
#' @examples 
#' result <- matrix.tau(c(1,2,3), c(2,3,4,5,6))
#' print(result)
#' 
#' @export
matrix.tau <- function(x, y) {
  
  # check arguments
  if (missing(x)) stop('Missing input vector.')
  if (missing(y)) y <- x
  
  # compute x.0, an arbitrary start point for vector x.
  # This improves accuracy if the offset is large
  # compared to the range of values
  x.0 <- min(x, y)
  x <- x - x.0
  y <- y - x.0
  
  # compute the time lag matrix
  tau <- outer(x, y, "-")
  
  # return to calling function
  return(as.matrix(tau))  
}

# -----------------------------------------------------------
# dcf
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

#' Compute the Discrete Correlation Function
#' 
#' \code{dcf} returns the Discrete Correlation Function estimates.
#'
#' @param tau (vector) list of lags at which to compute the CCF.
#' @inheritParams cross.correlate
#'
#' @return
#' A data frame containing columns:
#'  \item{tau}{the centre of the lag bins (vector)}
#'  \item{ccf}{the correlation coefficent in each lag bin}
#'  \item{n}{the number of data point pairs included in each lag bin}
#'
#' @section Notes:
#' Input are two time series (data frames with columns: \code{t, y, dy}
#' [optional]) Output is the correlation coefficient \code{ccf[i]} in different
#' lag bins \code{tau[i]}.
#'
#' In what follows we refer to the \code{t, y} values of time series 1
#' (\code{ts.1}) as \code{t.1, y.1}, and similarly for time series 2.
#' 
#' We first subtract the mean values from \code{y.1} and \code{y.2}. Then, 
#' within the \code{i}th lag bin (\code{tau[i] - dtau/2}, \code{tau[i] +
#' dtau/2}) we collect all pairs (\code{y.1, y.2}) of data for which \code{t.1 -
#' t.2} falls within the lag bin. Using these pairs of data we compute the sum
#' of their product and normalise it:
#' 
#' \code{cov[i] = (1/n) * sum_{j,k=1..n} (y.1[j] * y.2[k])}
#' 
#' Here, \code{j, k} are index arrays of length \code{n} that specify the points
#' of time series 1 and 2 (respectively) which pair-up within lag bin \code{i}. 
#' This gives a covariance. If \code{cov = TRUE} we keep these values. Otherwise
#' (default: \code{cov = FALSE}) we normalise by the product of the standard 
#' deviations of \code{y.1} and \code{y.2}.
#' 
#' The number of output lags \code{tau} is \code{2*lag.bins+1}, and the lags are
#' centred on zero. So they run from \code{-max.lag} to \code{+max.lag}. Any lag
#' bins containing fewer than min.pts pairs of points will have \code{ccf = NA}.
#'
#' If 'errors' are suppled for either or both time series (dx.1, dx.2) 
#' and \code{use.errors = TRUE} then the variance used 
#' in the denominator terms of correlation coefficient calculation will be
#' the 'excess variance', i.e. the total sample variance minus the 
#' mean square error. If dy.1 and/or dy.2 are a single number, this is 
#' assumed to be the same error for each data point. If \code{use.errors = FALSE}
#' (default) then the usual sample variance will be used.
#' 
#' @seealso \code{\link{cross.correlate}}, \code{\link{iccf}}
#' @examples 
#' ## Example using NGC 5548 data
#' res <- dcf(cont, hbeta, tau = seq(-200, 200, by = 5))
#' plot(res$tau, res$ccf, type = "l", col = "blue", lwd = 3, bty = "n")
#' grid()
#' 
#' @export
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
      warning(paste('CCF outside of (-1, +1): ', r.max))
  }
  
  # return final output
  return(data.frame(tau = tau, ccf = dcf, n = n.i))
  
}
