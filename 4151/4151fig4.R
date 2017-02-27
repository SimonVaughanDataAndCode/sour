band <- c("BAT","X4","X3","X2","X1","UVW2","UVM2","UVW1","U","B","V")
cols <- c("#e22bdf","blue","green","darkorange","red","black","#e22bdf","blue","green","darkorange","red")
ref <- "UVW2"
ns <- 100
hc <- F

xl <- c(-5.84,10.17)
if (ref=="UVW2") xl <- c(-10.17,5.84)
ccfa <- rep(0,12*349)
dim(ccfa) <- c(12,349)
tt <- 0

if (hc) postscript(paste("../",ref,".eps",sep=""),width=3,height=9)
par(mfrow=c(length(band),1),mgp=c(3,0.8,0),mar=c(0,4,0,2.5),
    oma=c(3,0.1,0.1,0.1))

# read in data and setup reference band
bb <- read.csv("~/R/sour/4151/4151.csv",as.is=T)
ts.2<- bb[bb$Filter==ref,1:3]
if ( tt > 0 ) ts.2$y <- ts.2$y - rmed(ts.2$t,ts.2$y,tt)

t0 <- Sys.time()
for (i in 1:length(band)) {
#for (i in c(1,5)) {
  ts.1 <- bb[bb$Filter==band[i],1:3]
  if ( tt > 0 ) ts.1$y <- ts.1$y - rmed(ts.1$t,ts.1$y,tt)
  
  iccf.result <- cross.correlate(ts.1, ts.2, method="iccf", nsim=ns, chatter=1, 
                                 peak.frac=0.8, dtau=0.1, max.lag=17.4)
  if(!all(is.finite(iccf.result$cent.dist))) {
    mask <- which(!is.finite((iccf.result$cent.dist)))
    cat(iccf.result$cent.dist[mask])
  }
  dtau <- diff(iccf.result$tau[1:2])
  brks <- (-300:300)*dtau
  hr <- hist(iccf.result$cent.dist, breaks=brks, main="",xaxt="n",col=cols[i],border=cols[i],
       yaxs="i",xaxs="i",xlim=xl*0.926,ylab=paste(band[i],"/",ref),las=1)

# compute the ICCF
#  lines(iccf.result$tau, iccf.result$ccf, col="black",lw=2)
  lines(iccf.result$tau, iccf.result$ccf*max(hr$counts), col="black",lw=2)
  
  abline(v=-100:100,col="grey")
  abline(v=-100:100,lty=2)
  abline(v=0,lwd=2)
  abline(h=max(hr$counts)/2,lty=3)
  axis(4,at=c(0,max(hr$counts)/2,max(hr$counts)),labels=c("","0.5","1"),las=1)
  box(lwd=2)
#  cat("Band ",i,": ",band[i],", r_max: ",mr,": ",band[i],", time ",Sys.time()-t0,fill=T,sep="")
  cat(i,band[i],max(iccf.result$ccf),quantile(iccf.result$cent.dist,probs=c(0.5,0.1,0.9)),
     min(iccf.result$tau),max(iccf.result$tau),length(iccf.result$tau),Sys.time()-t0,fill=T)
  
# write out CCF
  if(i==1) ccfa[1,] <- c("tau",iccf.result$tau)
  ccfa[i+1,] <- c(band[i],iccf.result$ccf)
}
write.csv(ccfa,file="ccfarray.csv",row.names=F)
axis(1)
axis(4,at=0,las=1)
mtext("Lag (days)",side=1,line=1.7,cex=0.7)
if (hc) graphics.off() 
