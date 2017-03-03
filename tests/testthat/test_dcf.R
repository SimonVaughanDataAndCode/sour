# tests for ICCF

context("DCF unit tests")

# check the difference matrix
  test_that("matrix_tau returns a matrix of the right form", {
    x <- seq(0,1,by=0.1)
    matr <- matrix_tau(x, x+0.5)
    expect_equal(mode(matr), "numeric")
    expect_equal(length(dim(matr)), 2)
  })

# Examples from Venables & Ripley
  require(graphics)
  tsf <- data.frame(t = time(fdeaths), y = fdeaths)
  tsm <- data.frame(t = time(mdeaths), y = mdeaths)
  tau <- seq(0, 1.5, by = 1/12)
  tau <- c(-rev(tau), tau[-1])

# compute CCF using raw ICCF method (unnormalised)
  ccf <- dcf(tsm, tsf, tau = tau)
  acf <- dcf(tsm, tsm, tau = tau)
  
  # compute CCF using DCF method 
  ccf <- dcf(tsm, tsf, tau = tau)
  acf <- dcf(tsm, tsm, tau = tau)
  
  test_that("dcf returns a list with tau and ccf columns", {
    expect_equal(mode(ccf), "list")
    expect_true(exists("tau", ccf))
    expect_true(exists("ccf", ccf))
  })
  
  # NOTE: ACF computed with DCF doesn't always give ACF(0)=1.
  test_that("DCF computed ACF[i] peaks at i=1.", {
    expect_equal(which.max(acf$ccf), 19)
  })
  
  