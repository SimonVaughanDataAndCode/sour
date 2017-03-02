#' NGC 5548 reverberation mapping data 1988-2001.
#'
#' Ground-based optical light curves of NGC 5548. Listed are the 5100A 
#' continuum fluxes. 
#'
#' @docType data
#'
#' @usage data(cont)
#'
#' @format A data frame with 1548 rows and 3 variables.
#'   \describe{
#'     \item{t}{The time (in MJD).}
#'     \item{y}{The flux.}
#'     \item{dy}{The (1-sigma) error on the flux.}
#'   }
#'
#' @keywords datasets
#'
#' @references B. Peterson et al., 2001, ApJ, v581, pp197, 204.
#' (\href{http://adsabs.harvard.edu/abs/2002ApJ...581..197P}{ADS})
#'
#' @source \href{http://www.astronomy.ohio-state.edu/~agnwatch/n5548/lcv/}{AGN Watch data}
#'
#' @examples
#' plot(cont$t, cont$y, bty = "n", type = "o", xlab = "time (day)", ylab = "flux")
#' grid()
"cont"
