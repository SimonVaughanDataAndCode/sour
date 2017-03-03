# tests for ICCF

context("ICCF unit tests")

# Examples from Venables & Ripley
  require(graphics)
  tsf <- data.frame(t = time(fdeaths), y = fdeaths)
  tsm <- data.frame(t = time(mdeaths), y = mdeaths)
  tau <- seq(0, 1.5, by = 1/12)
  tau <- c(-rev(tau), tau[-1])

# compute CCF using raw ICCF method (unnormalised)
  ccf <- iccf_core(tsm$y, tsm$y, tsf$t, tsf$y, tau = tau)
  acf <- iccf_core(tsm$t, tsm$y, tsm$t, tsm$y, tau = tau)
  
  test_that("iccf_core returns a list with an r column", {
    expect_equal(mode(ccf), "list")
    expect_true(exists("r", ccf))
  })
  
  test_that("ICCF computed ACF[i] peaks at i=1.", {
    expect_equal(which.max(acf$r), 19)
  })
  
  # compute CCF using full ICCF method (normalised)
  ccf <- iccf(tsm, tsf, tau = tau)
  acf <- iccf(tsm, tsm, tau = tau)
  
  test_that("iccf returns a list with tau and ccf columns", {
    expect_equal(mode(ccf), "list")
    expect_true(exists("tau", ccf))
    expect_true(exists("ccf", ccf))
  })
  
  test_that("DCF computed ACF[tau] peaks at tau=0.", {
    expect_equal(which.max(acf$ccf), 19)
    expect_equal(acf$ccf[19], 1)
  })
  
  