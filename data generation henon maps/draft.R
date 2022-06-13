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
henonmaps_data <- henon_map(b,d,r[[i]],timelength+10000)
x <- henonmaps_data$x[10001:timelength+10000] # discard the first 10000 samples
y <- henonmaps_data$y[10001:timelength+10000]


