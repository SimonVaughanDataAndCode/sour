# function to measure a running median of width w
# -----------------------------------------------------------
rmed <- function(t,y,w) {
  res <- y*0
  zz <- y*0
  for (ii in 1:length(res)) {
    tt <- t[ii]
    yy <- y[t>(tt-w/2) & t<(tt+w/2)]
    res[ii] <- mean(yy)
    zz[ii] <- length(yy)
  }
#  plot(t,y)
#  lines(t,res)
#  cat(zz,fill=T)
  return(res)
}