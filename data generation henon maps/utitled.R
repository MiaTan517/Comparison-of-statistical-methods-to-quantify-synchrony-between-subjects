library(IRISSeismic)
henon_map <- function(b, d, mu, n) {
  x <- numeric()
  y <- numeric()
  x[1]<- runif(1)
  x[2] <- runif(1)
  y[1]<- runif(1)
  y[2] <- runif(1)
  for (k in 2:n) {
    x[k+1] <- 1.4+b*x[k-1]-(x[k])^2
  }
  
  for (k in 2:n) {
    y[k+1] <- 1.4+d*y[k-1]-(mu*x[k]+(1-mu)*y[k])*y[k]
  }
  return (list(x=x,y=y))
}
c <- c();f <- c();d <- c();cc <- c();ccn <- c();fn <- c()
kullback <- j <- mi <- da <- jr <- c()
for (r in seq(0,1,0.01)) {
  henonmaps_data <- henon_map(0.3,0.3,r,10054+10000)
  x <- henonmaps_data$x[10001:(10054+10000)] # discard the first 10000 samples
  y <- henonmaps_data$y[10001:(10054+10000)]
  out <- cwt_signal(Sig=x,Cwt='MORLET',Min=4,Max=30,
                    FreqEch=200,Senuil=0,frqsmp=1,Aff=0,
                    Fc=0.8125)
  ascwt <- out$ascwt
  ccfs <- out$ccfs
  s <- out$scales
  
  #Normalization
  for (k in 1:nrow(ascwt)) {
    ascwt[k,] <- ascwt[k, ]/mean(ascwt[k, ])
  }
  ascwt1 <- ascwt/sum(ascwt) #divided by total sum
  
  
  
  out <- cwt_signal(Sig=y,Cwt='MORLET',Min=4,Max=30,
                    FreqEch=200,Senuil=0,frqsmp=1,Aff=0,
                    Fc=0.8125)
  
  ascwt <- out$ascwt
  ccfs <- out$ccfs
  s <- out$scales
  # Normalization
  for (k in 1:nrow(ascwt)) {
    ascwt[k, ] <- ascwt[k, ]/mean(ascwt[k, ])
  }
  ascwt2 <- ascwt/sum(ascwt)
  
  # part 3
  v1 <- seq(0.1, 0.9, 0.4)
  v2 <- seq(2, 10, 4)
  inftheo <- divergence_measures(ascwt1, ascwt2, v1, v2)

  kullback <- c(kullback,inftheo$K)
  j <- c(j,inftheo$J)
  mi <- c(mi,inftheo$MI)
  da <- c(da,inftheo$Da)
  jr <- c(jr,inftheo$JR)
  c=c(c,cor(x,y))
  f=c(f,filtdir(x,y))
  ts1 <- ts(x)
  ts2 <- ts(y)
  DF <- crossSpectrum(ts.union(ts1,ts2), spans=c(5,7))
  d=c(d,mean(DF$coh))
  #varofnoise1 <- Varnoise(x,((1-0.2)/0.2)) 
  #varofnoise2 <- Varnoise(y,((1-0.2)/0.2)) 
  #noise1 <- rnorm(length(x),0,sqrt(varofnoise1)) 
  #noise2 <- rnorm(length(x),0,sqrt(varofnoise2))
  #xnoise <- x + noise1 #add noise to signal
  #ynoise <- y + noise2
  #fn=c(fn,filtdir(xnoise,ynoise))
}

r=seq(0,1,0.01)
d
f[62:68]
c[59:64]
ccn
fn
min(c)
c;cc
f
which.min(c)
Varnoise <- function(signal, snr) varnoise <- mean(signal^2)/snr
varofnoise1 <- Varnoise(x,((1-0.2)/0.2)) 
varofnoise2 <- Varnoise(y,((1-0.2)/0.2)) 
noise1 <- rnorm(length(x),0,sqrt(varofnoise1)) 
noise2 <- rnorm(length(x),0,sqrt(varofnoise2))
xnoise <- x + noise1 #add noise to signal
ynoise <- y + noise2
c[11:30];f[11:30];d[11:30]
c[41:60];f[41:60];d[41:60]
plot(r[11:30],c[11:30],type='l',ylim=(c(-0.03,0.3)))
lines(r[11:30],f[11:30],col='red')
lines(r[11:30],d[11:30],col='orange')
lines(r[11:30],rep(1,20)*d[11])
lines(r[11:30],rep(1,20)*kullback[11])
plot(r[11:30],mi[11:30])

plot(r,c,type='l')
lines(r,f,col='red')
lines(r,d,col='orange')
lines(r,j)
lines(r,kullback,col='orange')
lines(r,jr,col='orange')
lines(r,d,col='orange')
lines(r[11:30],kullback[11:30])
kullback[41:60]
lines(r,r,col='blue')
lines(r,abs(c),col='purple')
lines(r,ccn,col='red')
plot(x,y)
fs <- 10
ts1 <- ts(x,frequency = fs)
ts2 <- ts(y, frequency = fs)
DF <- crossSpectrum(ts.union(ts1,ts2), spans=c(5,7))
mean(DF$coh)

y=0.1*x
cross_correlation(x,y,0,length((y)))
filtdir(a,b)
synchro <- function(x, m, tau, theiler, nn) {
  ndatos <-  nrow(x)
  xn <- scale(x)
  # ccf(x=x[,1], y=x[,2],lag.max = 0,type='correlation', plot=F)
  cross <- mean(xn[,1]*xn[,2]) # biased cross-correlation at lag 0
  out <- c()
  rxy <- ryx <- 0
  sxy = 0; syx = 0; hxy = 0; hyx = 0; nxy = 0; nyx = 0
  auxx <- auxy <- indexx<- indexy<- distx <- disty <-  c()
  for (i in 1:(ndatos-(m-1)*tau-1)) {
    for (k in 1:nn) { 
      # INICIALIZE AUX 
      auxx[k] <- 100000000
      indexx[k] <- 100000000
      auxy[k] <- 100000000
      indexy[k] <- 100000000
    }
    auxx[nn+1] <- 0
    auxy[nn+1] <- 0
    indexx[nn+1] <- 100000000
    indexy[nn+1] <- 100000000
    rrx = 0; rry = 0
    for (j in 1:(ndatos-(m-1)*tau-1)) {
      distx[j] <- 0
      disty[j] <- 0
      for (k in 0:(m-1)) {
        # DISTANCES
        distx[j] <- distx[j] + (x[i+k*tau, 1] - x[j+k*tau,1])^2
        disty[j] <- disty[j] + (x[i+k*tau, 2] - x[j+k*tau,2])^2
      }
      if ((abs(i-j) > theiler)) {
        if (distx[j]<auxx[1]) {
          flagx <- 0
          for (k in 1:(nn+1)) {
            if (distx[j]<auxx[k]) {
              auxx[k] <- auxx[k+1]
              indexx[k] <- indexx[k+1]
            } else  {
              auxx[k-1] <- distx[j]
              indexx[k-1] <- j
              flagx <- 1
            }
            if (flagx==1) break
          }
        }
        if (disty[j]<auxy[1]) {
          flagy <- 0
          for (k in 1:(nn+1)) {
            if (disty[j]<auxy[k]) {
              auxy[k] <- auxy[k+1]
              indexy[k] <- indexy[k+1]
            } else  {
              auxy[k-1] <- disty[j]
              indexy[k-1] <- j
              flagy <- 1
            }
            if (flagy==1) break
          }
        }
      }
      rrx <- rrx + distx[j]
      rry <- rry + disty[j]
      
    }
    rxx <- 0; ryy <- 0; rxy <- 0; ryx <- 0
    for (k in 1:nn) {
      rxx <- auxx[k] + rxx
      ryy <- auxy[k] + ryy
      rxy <- distx[indexy[k]] + rxy
      ryx <- disty[indexx[k]] + ryx
    }
    rxx <- rxx/nn
    ryy <- ryy/nn
    rxy <- rxy/nn
    ryx <- ryx/nn
    
    hxy <- hxy + log(rrx/((ndatos-(m-1)*tau-2)) / rxy)
    hyx <- hyx + log(rry/((ndatos-(m-1)*tau-2)) / ryx)
   
  }
  
  
 
  hxy = hxy/(ndatos-(m-1)*tau-1)
  hyx = hyx/(ndatos-(m-1)*tau-1)
 
  
  out["sxy"]=sxy
  out["syx"]=syx
  out["hxy"]=hxy
  out["hyx"]=hyx
  out["nxy"]=nxy
  out["nyx"]=nyx
  out["cross"]=cross
  
  return (out)
  
}
synchro(cbind(x,10,0,50,10)
synchro()