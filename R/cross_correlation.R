# To do:
#   - iccf.core: bring approx() interpolation outside main loop (speed)
#   - add more test data, unit tests
#   - use ACF(0) to normalise CCF(0) to ensure ACF(0)=1...?
# -----------------------------------------------------------


#' Estimate cross correlation of unevenly sampled time series.
#' 
#' \code{cross_correlate} returns cross correlation data for two times series.
#' 
#' Function for estimating the cross-correlation between two time series which
#' may be irregularly and/or non-simultaneously sampled. The CCF is computed
#' using one of two methods: (1) the Discrete Correlation Function (DCF; Edelson
#' & Krolik 1988) or (2) the Interpolated Cross Correlation Function (ICCF; 
#' Gaskell & Sparke 1986). You can also produce estimates of uncertainty on the
#' CCF, its peak and centroid using the Flux Randomisation and Random Subsample
#' Selection (FR/RSS) method of Peterson et al. (1998).
#' 
#' @param ts.1,ts.2    (array or dataframe) data for time series 1 and 2.
#' @param method       (string) use \code{"dcf"} or \code{"iccf"} (default).
#' @param max.lag      (float) maximum lag at which to compute the CCF.
#' @param min.pts      (integer) each DCF bin must contain at least \code{min.pts} correlation coefficients.
#' @param dtau         (float) spacing of the time delays (\code{tau}) which which CCF is estimated.
#' @param local.est    (logical) use 'local' (not 'global') means and variances?
#' @param zero.clip    (logical) remove pairs of points with exactly zero lag?
#' @param use.errors   (logical) if \code{TRUE} then subtract mean square error from variances.
#' @param one.way      (logical) (ICCF only) if TRUE then only interpolar time series 2.
#' @param cov          (logical) if \code{TRUE} then compute covariance, not correlation coefficient.
#' @param prob         (logical) probability level to use for confidence intervals
#' @param nsim         (integer) number of FR/RSS simulations to run 
#' @param peak.frac    (float) only include CCF points above \code{peak.frac}*max(ccf) in centroid calculation.   
#' @param chatter      (integer) set the level of feedback.
#' @param plot         (logical) if \code{TRUE} then a plot of the ccf vs. tau is produced.
#' @param ...          (other) any other plot function parameters.
#' 
#' @return 
#'  A list with components
#'  \item{tau}{(array) A one dimensional array containing the lags at which the CCF is estimated.}
#'  \item{ccf}{(array) An array with the same dimensions as lag containing the estimated CCF.}
#'  \item{lower}{(array) Lower limit of CCF (see Notes).}
#'  \item{upper}{(array) Upper limit of CCF (see Notes).}
#'  \item{peak.dist}{(array) A array of length \code{nsim} containing the CCF peaks from the simulations.}
#'  \item{cent.dist}{(array) A array of length \code{nsim} containing the CCF centroids from the simulations.}
#'  \item{method}{(string) which method was used? \code{"iccf"} or \code{"dcf"}}
#'  The value of \code{ccf[k]} returns the estimated correlation between
#'  \code{ts.1$y}(t+tau) and \code{ts.2$y}(t) where tau = \code{tau[k]}. A
#'  strong peak at negative lags indicates the \code{ts.1} leads \code{ts.2}.
#'  
#' @section Notes:
#' If only one time series is given as input then the Auto-Correlation Function
#' (ACF) is computed.
#' 
#' Input data frames: note that the input data \code{ts.1} and \code{ts.2} are 
#' not traditional \code{R} time series objects. Such objects are only suitable 
#' for regularly sampled data. These CCF functions are designed to work with 
#' data of arbitrary sampling, we therefore need to explicitly list times and 
#' values. The input objects are therefore data.frames with at least two columns
#' which much be called \code{t} (time) and \code{y} (value). An error of the 
#' value may be provided by a \code{dy} column. Any other columns are ignored.
#' 
#' Local vs. global estimation: If \code{local.est = FALSE} (default) then the
#' correlation coefficient is computed sing the 'global' mean and variance of
#' each time series. If \code{local.est = TRUE} then the correlation coefficient
#' is computed using the 'local' mean and variance. For each lag, the mean to be
#' subtrated and the varaince to be divided are computed using only data points
#' contributing to that lag.
#' 
#' Simulations: Performs "flux randomisation" and "random sample selection" of 
#' an input time series, following Peterson et al. (2004, ApJ, 613:682-699).
#' Given an input data series \code{(t, y, dy)} of length \code{N} we sample 
#' \code{N} points with replacement. Duplicated points are ignored, so the 
#' ouptut is usually shorter than the input. So far this is a basic bootstrap 
#' procedure.
#' 
#' If error bars are provided: when a point is selected \code{m} times, we
#' decrease the error by \code{1/sqrt(m)}. See Appendix A of Peterson et al. And
#' after resampling in time, we then add a random Gaussian deviate to each
#' remaining data point, with std.dev equal to its error bar. In this way both
#' the times and values are randomised. If errors bars are not provided, this is
#' a simple bootstrap.
#'
#' Peak and centroid: from the simulations we record the CCF peak and its
#' centroid. The centroid is the mean of \code{tau*ccf/sum(ccf)} including all
#' points for which \code{ccf} is higher than \code{peak.frac} of
#' \code{max(ccf)}.
#'  
#' Upper/lower limits and distributions: If no simulations are used (\code{nsim
#' = 0}) then upper/lower confidence limits on the CCF are estimated using the
#' method of Barlett (1955) based on the two ACFs. If simulations are used, the
#' confidence limits are based on the simulations.
#' 
#' @seealso \code{\link[stats]{ccf}}, \code{\link{fr_rss}}
#' 
#' @examples
#'  ## Example using NGC 5548 data
#'  result <- cross_correlate(cont, hbeta, method = "iccf", dtau = 1, max.lag = 550)
#'  plot(result$tau, result$ccf, type = "l", bty = "n", xlab = "time delay", ylab = "CCF")
#'  grid()
#'  
#'  ## or using the DCF method
#'  result <- cross_correlate(cont, hbeta, method = "dcf", dtau = 5, max.lag = 350)
#'  lines(result$tau, result$ccf, col = "red")
#'
#'  ## Examples from Venables & Ripley
#'  require(graphics)
#'  tsf <- data.frame(t = time(fdeaths), y = fdeaths)
#'  tsm <- data.frame(t = time(mdeaths), y = mdeaths)
#' 
#'  ## compute CCF using ICCF method
#'  result <- cross_correlate(tsm, tsf, plot = TRUE, method = "iccf")
#' 
#'  ## compute CCF using standard method (stats package) and compare
#'  result.st <- ccf(mdeaths, fdeaths, plot = FALSE)
#'  lines(result.st$lag, result.st$acf, col="red", lwd = 3)
#' 
#' @export
cross_correlate <- function(ts.1, ts.2,
                            method = "iccf",
                            max.lag = NULL, 
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
    sims <- ccf_errors(ts.1, ts.2, tau, nsim = nsim,
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
      acf.2 <- dcf(ts.2, ts.2, tau,
                   local.est = local.est, 
                   min.pts = min.pts, 
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
# fr_rss
# History:
#  21/03/16 - First working version
#
# (c) Simon Vaughan, University of Leicester
# -----------------------------------------------------------

#' Perform flux randomisation/random subset section on input data.
#' 
#' \code{fr_rss} returns a randomise version of an input data array.
#' 
#' Performs "flux randomisation" and "random sample selection"
#' of an input time series, following Peterson et al. (2004, ApJ, v613, pp682-699).
#' This is essentially a bootstrap for a data vector.
#' 
#' @param dat (data frame) containing columns \code{t, y} and (optionally) 
#'            \code{dy}.
#'
#' @return
#'  A data frame containing columns
#'  \item{t}{time bins for randomised data}
#'  \item{y}{values for randomised data}
#'  \item{dy}{errors for randomised data}
#'
#' @section Notes:
#' Given an input data series (\code{t, y, dy}) of length \code{N} we sample 
#' \code{N} points with replacement. Duplicated points are ignored, so the
#' ouptut is usually shorter than the input. So far this is a basic bootstrap
#' procedure.
#' 
#' If error bars are provided: when a point is selected \code{m} times, we
#' decrease the error, scaling by \code{1/sqrt(m)}. See Appendix A of Peterson
#' et al. After resampling, we then add a random Gaussian deviate to each
#' remaining data point, with std.dev equal to its (new) error bar. If errors
#' bars are not provided, this is a simple bootstrap (no randomisation of
#' \code{y}).
#' 
#' @seealso \code{\link{cross_correlate}}, \code{\link{ccf_errors}} 
#' 
#' @examples 
#'  ## Example using the NGC 5548 data
#'  plot(cont$t, cont$y, type="l", bty = "n", xlim = c(50500, 52000))
#'  rcont <- fr_rss(cont)
#'  lines(rcont$t, rcont$y, col = "red")
#' 
#'  ## Examples from Venables & Ripley
#'  require(graphics)
#'  plot(fdeaths, bty = "n")
#'  tsf <- data.frame(t = time(fdeaths), y = fdeaths)
#'  rtsf <- fr_rss(tsf)
#'  lines(rtsf$t, rtsf$y, col="red", type="o")
#'
#' @export
fr_rss <- function(dat) {
  
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
# ccf_errors
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

#' Use random simulations to estimate undertainty on CCF estimates.
#' 
#' \code{ccf_errors} returns information on the undertainty on CCF estimates.
#'
#' Computes errors on the CCF estimates using "flux randomisation" 
#' and "random subset sampling" FR/RSS using the \code{fr_rss} function.
#'
#' @param tau (array) list of lags at which the CCF is to be evaluated.
#' @param acf.flag (TRUE)logical) \code{TRUE} when computing ACF, and \code{ts.2 = ts.1}
#' @inheritParams cross_correlate
#'
#' @return 
#' The output is a list containing two data frames: \code{lags} and \code{dists}.
#'    \item{lags}{a data frame with four columns}
#'    \item{tau}{time lags}
#'    \item{dcf}{the DCF values for the input data}
#'    \item{lower}{the lower limit of the confidence interval}
#'    \item{upper}{the upper limit of the confidence interval}
#'    \item{dists}{a data frame with two columns}
#'    \item{peak.lag}{the peak values from nsim simulations}
#'    \item{cent.lag}{the centroid values from nsim simulations}
#'
#' @section Notes:
#' For each randomised pair of light curves we compute the CCF. We record the
#' CCF, the lag at the peak, and the centroid lag (including only points higher
#' than \code{peak.frac * max(ccf)}). Using \code{nsim} simulations we compute
#' the \code{(1-p)*100\%} confidence intervals on the CCF values, and the
#' distribution of the peaks and centroids.
#'
#' @seealso \code{\link{cross_correlate}}, \code{\link{fr_rss}} 
#'
#' @export
ccf_errors <- function(ts.1, ts.2, 
                       tau = NULL,
                       min.pts = 5,
                       local.est = FALSE,
                       cov = FALSE,
                       prob = 0.1, 
                       nsim = 250,
                       peak.frac = 0.8,
                       zero.clip = NULL,
                       one.way = FALSE,
                       method = "iccf",
                       use.errors = FALSE,
                       acf.flag = FALSE,
                       chatter = 0) {

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
    ts1.sim <- fr_rss(ts.1)
    if (acf.flag == FALSE) {
      ts2.sim <- fr_rss(ts.2)
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
