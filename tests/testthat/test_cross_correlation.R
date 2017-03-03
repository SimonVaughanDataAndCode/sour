# tests for ICCF

context("CROSS_CORRELATION unit tests")

# Examples from Venables & Ripley
  require(graphics)
  tsf <- data.frame(t = time(fdeaths), y = fdeaths)
  tsm <- data.frame(t = time(mdeaths), y = mdeaths)
  tau <- seq(0, 1.5, by = 1/12)
  tau <- c(-rev(tau), tau[-1])

# compute CCF using  ICCF method 
  ccf <- cross_correlate(tsm, tsf, dtau = 1/12)
  acf <- cross_correlate(tsm, tsm, dtau = 1/12)
  
  test_that("cross_correlate returns a list with tau and ccf columns", {
    expect_equal(mode(ccf), "list")
    expect_true(exists("tau", ccf))
    expect_true(exists("ccf", ccf))
  })
  
  # NOTE: ACF computed with DCF doesn't always give ACF(0)=1.
  test_that("CCF computed ACF[i] peaks at tau = 0.", {
    expect_equal(which.max(acf$ccf), 19)
  })
  
  