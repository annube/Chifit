

sigma <- function(x) sqrt(1-4/x)
g <- function(x) sigma(x) * log( (sigma(x) - 1) / ( sigma(x) + 1 ) ) + 2

integrand <- function(y,x,r,k){
  print(paste("eval ", length(y)))
  re = k %% 2 == 0
  if(re) y^k * exp(-x*sqrt(1+y^2) ) * 2* Re( g(2+2i*r*y)  )
  else y^k * exp(-x*sqrt(1+y^2) ) *  2*Im( g(2+2i*r*y)  )
}



k=0
r=0.82
x0=4

plot( function(x) integrand(x,x0,r,k) , xlim=c(0,5) )


int.res <- integrate(function(x) integrand(x,x0,r,k) ,
                     lower=0,upper=Inf,
                     rel.tol=1.e-6,
                     abs.tol=0
                     )



## use gamma distribution for the mapping
## this is a good crosscheck but needs much more subdivisions in the numeric integrator
## and as we also need to evaluate qgamma and dgamma I think the performance
## is much worse

integrand.gamma <- function(z,x,r,k) {
  y <- qgamma(z,shape=k+1,rate=x)
  integrand(y,x,r,k)/dgamma(y,shape=k+1,rate=x)
}

plot( function(x) integrand.gamma(x,x0,r,k) , xlim=c(0,1) ,n=500)

int.res.gamma <- integrate(function(x) (integrand.gamma(x,x0,r,k)) ,
                           lower=0,upper=1,
                           rel.tol=1.e-3,
                           abs.tol=0
                           )
