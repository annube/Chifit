

sigma <- function(x) sqrt(1-4/x)
g <- function(x) sigma(x) * log( (sigma(x) - 1) / ( sigma(x) + 1 ) ) + 2

integrand <- function(y,x,r,k){
  y^k * exp(-x*sqrt(1+y^2) ) * ( g(2+2i*r*y) + (-1)^k * g(2-2i*r*y) )
}



k=0
r=1
x0=10

plot( function(x) integrand(x,x0,r,k) , xlim=c(0,5) )


int.res <- integrate(function(x) Re(integrand(x,x0,r,k)) , lower=0,upper=Inf)


integrand.gamma <- function(z,x,r,k) {
  print("eval")
  y <- qgamma(z,shape=k+1,rate=x)
  integrand(y,x,r,k)/dgamma(y,shape=k+1,rate=x)
}

plot( function(x) integrand.gamma(x,x0,r,k) , xlim=c(0,1) ,n=500)

int.res.gamma <- integrate(function(x) Re(integrand.gamma(x,x0,r,k)) , lower=0,upper=1)
