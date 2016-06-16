# Sour
Functions for cross correlation of unevenly sampled time series.

These are pure R functions for estimating cross correlation functions (CCF) based to two time series that may be unevenly and asynchronously sampled time series. Options are the Discrete Correlation Function (DCF; Edelson & Krolik 1988) or the Interpolated Cross Correlation Function (ICCF; Gaskell & Sparke 1986). 
You can also produce estimates of uncertainty on the CCF, its peak and centroid using the Flux Randomisation and Random Subsample Selection (FR/RSS) method of Peterson et al. (1998). 

## Setting up

Just load the functions

```R
 source("cross_correlation.R")
```

and you should be ready to use the main function cross.correlate(...)

For a worked example: see the tests/test_CrossCorrelation.R script

## Example

![example](figures/time_series.png)

Example of two time series to compare. Each series is an array or data frame with columns t, y (and optionally) dy. In this example, we have ts.1 and ts.2 as the two data frames. 

```R
  ts.1 <- read.table("data/time_series_1.txt", header=TRUE)
  ts.2 <- read.table("data/time_series_2.txt", header=TRUE)

  # insert some dodgy (missing) data
  ts.1[20:50,2] <- NA
  ts.2[40:70,2] <- NA
```

Then we compute the DCF with

```R
  dcf.out <- cross.correlate(ts.1, ts.2, local.est = TRUE, 
                             dtau = 1, nsim = 2000, max.lag = 70)
```

Here we chose the width for the lag bins (1.0), the maximum lag to examine (-70.0 to +70.0), and use 2,000 simulations to estimate the centroid distribution. Setting local.est = TRUE means that the mean and variances are computed using only pairs of data contributing to a given lag bin.

![example](figures/ccf.png)

The resulting CCF computed with the DCF. The is plotted using, e.g.

```R
  plot(dcf.out$tau-dtau/2, dcf.out$ccf, type = "s", lwd = 2, 
       bty = "n", ylim = c(-1, 1), xlab = "lag", ylab = "CCF", 
       col = "black", main = "DCF test - true lag is 10")
  grid(col = "darkgrey")
```

![example](figures/ccf_centroid_distribution.png)

The distribution of centroids of the CCF from 2,000 simulations. Plotted using e.g.

```R
  hist(dcf.out$cent.dist, breaks = 50, col = "steelblue1", 
       border = NA, main = "CCF centroid distribution", prob = TRUE)
  lines(density(dcf.out$cent.dist, n = 256, na.rm = TRUE), lwd = 2)
  cat('-- mean lag', mean(dcf.out$cent.dist), fill = TRUE)
```

## References

For more info on the methods see:

R. Edelson & J. Krolik (1988; ApJ) http://adsabs.harvard.edu/abs/1988ApJ...333..646E

C. M. Gaskell & L. S. Sparke (1986; ApJ) http://adsabs.harvard.edu/abs/1986ApJ...305..175G

B. Peterson et al. (1998; PASJ) http://adsabs.harvard.edu/abs/1998PASP..110..660P

