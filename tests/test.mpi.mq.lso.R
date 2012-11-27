library(chifit)



## 
afpi.mq <- function(x,par,aargs){

  R = numeric(aargs$numLs)
  R[1] = 1
  if(aargs$numLs > 1 )
    R[2:aargs$numLs] = par[ 4 + 1:(aargs$numLs-1) ]
  
  print(R)
  
  aB0 <- par[1] * R[aargs$lsIndex]
  af0 <- par[2] * R[aargs$lsIndex]
  aLambda4 <- par[3] * R[aargs$lsIndex]
  Cf <- par[4] * R[aargs$lsIndex]^2
  chi_mu = 2*aB0*x

  xi.ll <-  chi_mu / ( 4 * pi * af0 )^2
  
  return(  af0 * ( 1 -
                     2 * xi.ll *
                     log(chi_mu / (aLambda4)^2 ) + Cf
                 ) * ( 1 -2 *xi.ll* g.tilde( sqrt( chi_mu ) * aargs$L ,partitions ) )
         )
  
}


x = c(
  0.01,
  0.005,
  0.008)

par = c(
  1.1,
  0.05,
  0.67,
  0.0,
  0.9
       )

aargs = list( numLs = 2,
  lsIndex = c(1,1,2),
  L = c(24,25,24),
  ZP = c(1,1,1)
  )


chi_mu = 2 * par[1] * x[3]

xi.ll = chi_mu / ( 4 * pi * par[2] ) ^ 2
