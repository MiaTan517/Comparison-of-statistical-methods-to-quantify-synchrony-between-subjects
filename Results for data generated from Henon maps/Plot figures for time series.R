library(ggplot2)
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
henonmaps_data <- henon_map(0.3,0.3,0.6,(1006+10000))
x <- henonmaps_data$x[10001:(1006+10000)] # discard the first 10000 samples
y <- henonmaps_data$y[10001:(1006+10000)]

ts1 <- data.frame(cbind(1:length(x),x))
ts2 <- data.frame(cbind(1:length(y),y))
ggplot(data = ts1[1:200,],mapping = aes(x=V1))+
  geom_line(aes(y=x))+scale_x_continuous(breaks = NULL)+xlab('Time length')+
  scale_y_continuous(limits = c(-2,2))+ylab('Amplitude of x')
ggplot(data = ts2[1:200,],mapping = aes(x=V1))+
  geom_line(aes(y=y))+scale_x_continuous(breaks = NULL)+xlab('Time length')+
  scale_y_continuous(limits = c(-2,2))+ylab('Amplitude of y')
##################################################################
library(BAMBI) # the bivariate von Mises sine distribution
n <- 161 #161 10054 17 1006
x<- c();y<-c()
a <- (1:n)*pi # generated an unwrapped and perfectly regular phase series 
#generated n independent samples from a von Mises random distribution
Varnoise <- function(signal, snr) varnoise <- mean(signal^2)/snr

a[1] <- 0
r=0.6

for (lambda in   2*log( (1+r) / (1-r) )) {
  print(lambda)
  b <- rvmsin(n,kappa1 = 2,kappa2 =2, kappa3 = lambda
              , mu1=pi, mu2=pi) #correlation=0.8, Kappa1,Kappa2= 8 
  b <- b-pi #adjust the range from [0,2*pi] to [-pi,pi]
  # add a and b together 
  s1 <- a+b[,1] #  
  s2 <- a+b[,2] #
  #We still need spline interpolation and take the sine of the data to get the simulated data
  f <- 10 # 
  omega <- 2*pi*f # 1/omega is the interval for spline interpolation
  #spline interpolation
  s1inter <- spline(1:length(s1), s1, xout = seq(1, length(s1), by = 1/omega))$y
  s2inter <- spline(1:length(s2), s2, xout = (seq(1, length(s2), by = 1/omega)))$y
  #take the sine
  x <- sin(s1inter); y <- sin(s2inter)
  varofnoise1 <- Varnoise(x,((1-0.5)/0.5)) 
  varofnoise2 <- Varnoise(y,((1-0.5)/0.5)) 
  noise1 <- rnorm(length(x),0,sqrt(varofnoise1)) 
  noise2 <- rnorm(length(x),0,sqrt(varofnoise2))
  xnoise <- x + noise1 #add noise to signal
  ynoise <- y + noise2
}



ts3 <- data.frame(cbind(1:length(x),x))
ts4 <- data.frame(cbind(1:length(y),y))
ggplot(data = ts3[1:2000,],mapping = aes(x=V1))+
  geom_line(aes(y=x))+scale_x_continuous(breaks = NULL)+xlab('Time length')+
  scale_y_continuous(limits = c(-2,2))+ylab('Amplitude of x')
ggplot(data = ts4[1:2000,],mapping = aes(x=V1))+
  geom_line(aes(y=y))+scale_x_continuous(breaks = NULL)+xlab('Time length')+
  scale_y_continuous(limits = c(-2,2))+ylab('Amplitude of y')


ts5 <- data.frame(cbind(1:length(xnoise),xnoise))
ts6 <- data.frame(cbind(1:length(ynoise),ynoise))
ggplot(data = ts5[1:2000],mapping = aes(x=V1))+
  geom_line(aes(y=xnoise))+scale_x_continuous(breaks = NULL)+xlab('Time length')+
  scale_y_continuous(limits = c(-2,2))+ylab('Amplitude of x')
ggplot(data = ts6[1:2000,],mapping = aes(x=V1))+
  geom_line(aes(y=ynoise))+scale_x_continuous(breaks = NULL)+xlab('Time length')+
  scale_y_continuous(limits = c(-2,2))+ylab('Amplitude of y')

