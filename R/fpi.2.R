# library(elliptic,pos="package:base")


fpi.2.calc.cache <- function(x,par,aargs) {
  
  nens <- aargs$num.ens
  pionrange <- 1:nens
  
  af0 <- par[aargs$ls.index[pionrange]+1]
  C= 4 * pi * af0 * aargs$Loa
  xi.ll.FV <- (x[pionrange]/(4*pi*af0[pionrange]))^2

  xi.ll <- numeric(nens)
  for( i in 1:nens ) {
    xi.ll[i] = determine.xi.IV(xi.ll.FV[i] , C[i])
  }
  
  xi.ll.times.g.tilde = xi.ll*g.tilde(sqrt(xi.ll)*C,partitions)

  cache=list( xi.ll = xi.ll ,
    xi.ll.times.g.tilde = xi.ll.times.g.tilde 
    )

  return(cache)

}


dfpi.2.calc.cache <- function(x,par,aargs) {
  
  nens <- aargs$num.ens
  pionrange <- 1:nens
  
  af0 <- par[aargs$ls.index[pionrange]+1]
  C= 4 * pi * af0 * aargs$Loa
  xi.ll.FV <- (x[pionrange]/(4*pi*af0[pionrange]))^2

  xi.ll <- numeric(nens)
  for( i in 1:nens ) {
    xi.ll[i] = determine.xi.IV(xi.ll.FV[i] , C[i])
  }
  
  xi.ll.times.g.tilde = xi.ll*g.tilde(sqrt(xi.ll)*C,partitions)
  xi.ll.times.g.tilde.deri = xi.ll * g.tilde.deri(sqrt(xi.ll)*C,partitions)

  dxi.ll.dC <- numeric(length(pionrange))
  for( i in pionrange ) {
    dxi.ll.dC[i] = ddetermine.xi.IV.dC( xi.ll.FV[i] , xi.ll[i], C[i],
               cache = list(
                 xi.ll.times.g.tilde = xi.ll.times.g.tilde[i],
                 xi.ll.times.g.tilde.deri=xi.ll.times.g.tilde.deri[i]
                 )
               )
  }

  
  cache=list(
    xi.ll = xi.ll ,
    xi.ll.times.g.tilde = xi.ll.times.g.tilde,
    xi.ll.times.g.tilde.deri = xi.ll.times.g.tilde.deri
    ,
    dxi.ll.dC = dxi.ll.dC
    
    )

  return(cache)

}

## chi-PT formula for the pion decay constant
fpi.2 <- function(x,par,aargs,cache=NULL){

  
  nens <- aargs$num.ens
  pionrange <- 1:nens
  
  b <- par[1]
  af0 <- par[aargs$ls.index[pionrange]+1]
  C= 4 * pi * af0 * aargs$Loa
  xi.ll.FV <- (x[pionrange]/(4*pi*af0[pionrange]))^2


  ## receycle precomputed xi.ll
  if( is.null(cache) ) {
    
    xi.ll <- numeric(nens)
    for( i in 1:nens ) {
      xi.ll[i] = determine.xi.IV(xi.ll.FV[i] , C[i])
    }

    xi.ll.times.g.tilde = xi.ll*g.tilde(sqrt(xi.ll)*C,partitions)
  } else {
   
    xi.ll = cache$xi.ll
    xi.ll.times.g.tilde = cache$xi.ll.times.g.tilde

  }

  Factor.1 = af0
  Factor.2 = ( 1-2*xi.ll*log(xi.ll) + xi.ll  *b )
  Factor.3 = (1-2 * xi.ll.times.g.tilde )
  
  return(
         1
         * Factor.1
         * Factor.2
         * Factor.3
         )
}


## chi-PT formula for the pion decay constant
dfpi.2 <- function(x,par,aargs,cache=NULL){

  
  nens <- aargs$num.ens
  pionrange <- 1:nens
  
  b <- par[1]
  af0 <- par[aargs$ls.index[pionrange]+1]
  C= 4 * pi * af0 * aargs$Loa
  xi.ll.FV <- (x[pionrange]/(4*pi*af0[pionrange]))^2
  

  ## receycle precomputed xi.ll
  if( is.null(cache) ) {
    
    xi.ll <- numeric(nens)
    for( i in 1:nens ) {
      xi.ll[i] = determine.xi.IV(xi.ll.FV[i] , C[i])
    }

    xi.ll.times.g.tilde = xi.ll*g.tilde(sqrt(xi.ll)*C,partitions)
    xi.ll.times.g.tilde.deri = xi.ll * g.tilde.deri(sqrt(xi.ll)*C,partitions)
    
    dxi.ll.dC <- numeric(length(pionrange))
    for( i in 1:nens ) {
      dxi.ll.dC[i] = ddetermine.xi.IV.dC( xi.ll.FV[i] , xi.ll[i], C[i]
                  ,
                  cache = list(
                    xi.ll.times.g.tilde = xi.ll.times.g.tilde[i],
                    xi.ll.times.g.tilde.deri=xi.ll.times.g.tilde.deri[i]
                    )
                 )
    }

  } else {
    xi.ll = cache$xi.ll
    xi.ll.times.g.tilde = cache$xi.ll.times.g.tilde
    xi.ll.times.g.tilde.deri = cache$xi.ll.times.g.tilde.deri
    dxi.ll.dC = cache$dxi.ll.dC

  }


  Factor.1 = af0
  Factor.2 = ( 1-2*xi.ll*log(xi.ll) + xi.ll  *b )
  Factor.3 = (1-2 * xi.ll.times.g.tilde )

  dFactor.1.daf0 = 1
  dFactor.2.daf0 = ( - 2 * (log(xi.ll) + 1 ) + b ) * dxi.ll.dC * C/af0
  dFactor.3.daf0 = - 2 * C / af0 * (
                                    dxi.ll.dC * (
                                                 xi.ll.times.g.tilde/xi.ll
                                                 + C * 0.5 * xi.ll.times.g.tilde.deri/sqrt(xi.ll)
                                                 )
                                    + xi.ll.times.g.tilde.deri * sqrt(xi.ll)
                                    )

  dFactor.2.db = xi.ll

  
  F12 = Factor.1 * Factor.2
  F23 = Factor.2 * Factor.3
  F13 = Factor.1 * Factor.3

  dfpi.db = Factor.1 * dFactor.2.db * Factor.3
  
  dfpi.daf0.raw = dFactor.1.daf0 * F23 + dFactor.2.daf0 * F13 + dFactor.3.daf0 * F12

  dfpi.daf0 = array(0, dim = c(nens, aargs$num.ls ))

  for( i in 1:nens )
    dfpi.daf0[i,aargs$ls.index[i]] = dfpi.daf0.raw[i]


  
  return(
         cbind(dfpi.db, dfpi.daf0)
         )
}

f.pi.exp <- 130.7
m.pi.exp <- 135.0

fpi.2.determine.a <- function(b,af0){
  C <- ( m.pi.exp / ( 4 * pi * af0 ) )^2
  f <- function(a) a * f.pi.exp - af0 * ( 1 - 2 * a^2 * C * log( a^2 * C ) + b * a^2 * C )
  df <- function(a) f.pi.exp - af0 * (  - 4 * a * C * ( log( a^2 * C ) + 1 ) + 2 * a * C * b )
  
  ##nr.res <- newton_raphson(0.00038,f, df ,100,tol=1.e-10)
  nr.res <- uniroot(f, c(0.0001,0.0005) ,tol=1.e-10)
##  print(paste(" number of newton iterations : " , nr.res$iter ))
    return(
           nr.res$root
           )

  
}


fK.2 <- function(x,par,aargs,cache=NULL){

  nens <- aargs$num.ens
  pionrange <- 1:nens
  mssrange=(1:aargs$num.mss)+nens
  
  b <- par[1]
  af0 <- par[aargs$ls.index[pionrange]+1]
  af0.i <- af0[mss.ensemble.index]

  f0Koaf0 = par[(1+aargs$num.ls)+1]
  fmKoaf0 = par[(1+aargs$num.ls)+2]

  b0K = par[(1+aargs$num.ls)+3]
  bmK = par[(1+aargs$num.ls)+4]



  C= 4*pi * af0 *aargs$Loa
  xi.ll.FV <- (x[pionrange]/(4*pi*af0))^2

  
  if( is.null(cache) ) {
    
    xi.ll <- numeric(length(pionrange))
    for( i in 1:aargs$num.ens ) {
      xi.ll[i] = determine.xi.IV(  xi.ll.FV[i], C[i])
    }

    xi.ll.times.g.tilde = xi.ll*g.tilde(sqrt(xi.ll)*C,partitions)
  } else {
    xi.ll = cache$xi.ll
    xi.ll.times.g.tilde = cache$xi.ll.times.g.tilde

  }

  
  xi.ss <- (x[mssrange])^2/(4*pi*af0.i)^2


  Factor.1 = af0.i * ( f0Koaf0 + fmKoaf0 * xi.ss )
  Factor.2 = (
              1
              - 0.75*(xi.ll*log(xi.ll))[mss.ensemble.index]
              + (b0K+bmK*xi.ss)*xi.ll[mss.ensemble.index]
              )
  Factor.3 = (1-3/4* xi.ll.times.g.tilde )[mss.ensemble.index]
  
##print(length(Factor.2))
  return(
##         xi.ll[mss.ensemble.index]
         1
         * Factor.1
         * Factor.2
         * Factor.3
         )

  
}

dfK.2 <- function(x,par,aargs,cache=NULL){

  nens <- aargs$num.ens
  pionrange <- 1:nens
  mssrange=(1:aargs$num.mss)+nens
  
  b <- par[1]
  af0 <- par[aargs$ls.index[pionrange]+1]
  af0.i = af0[mss.ensemble.index]

  f0Koaf0 = par[(1+aargs$num.ls)+1]
  fmKoaf0 = par[(1+aargs$num.ls)+2]

  b0K = par[(1+aargs$num.ls)+3]
  bmK = par[(1+aargs$num.ls)+4]



  xi.ll.FV <- (x[pionrange]/(4*pi*af0))^2
  C= 4*pi*af0 * aargs$Loa
  C.i = C[mss.ensemble.index]
  
  if( is.null(cache) ) {
    
    xi.ll <- numeric(length(pionrange))
    for( i in 1:aargs$num.ens ) {
      xi.ll[i] = determine.xi.IV( xi.ll.FV[i] , C[i])
    }

    xi.ll.times.g.tilde = xi.ll*g.tilde(sqrt(xi.ll)*C,partitions)
    xi.ll.times.g.tilde.deri = xi.ll * g.tilde.deri(sqrt(xi.ll)*C,partitions)

    dxi.ll.dC <- numeric(length(pionrange))
    for( i in 1:aargs$num.ens ) {
      dxi.ll.dC[i] = ddetermine.xi.IV.dC( xi.ll.FV[i] , xi.ll[i], C[i]
                 ,
                 cache = list(
                   xi.ll.times.g.tilde = xi.ll.times.g.tilde[i],
                   xi.ll.times.g.tilde.deri=xi.ll.times.g.tilde.deri[i]
                   )
                 
                 )
    }
    

  } else {
    xi.ll = cache$xi.ll
    xi.ll.times.g.tilde = cache$xi.ll.times.g.tilde
    xi.ll.times.g.tilde.deri = cache$xi.ll.times.g.tilde.deri
    dxi.ll.dC = cache$dxi.ll.dC
    
  }



  
  xi.ll.i <- xi.ll[mss.ensemble.index]
  xi.ss <- (x[mssrange])^2/(4*pi*af0.i)^2


  Factor.1 = af0.i * ( f0Koaf0 + fmKoaf0 * xi.ss )
  Factor.2 = (
              1
              - 0.75*(xi.ll*log(xi.ll))[mss.ensemble.index]
              + (b0K+bmK*xi.ss)*xi.ll[mss.ensemble.index]
              )
  Factor.3 = (1-3/4* xi.ll.times.g.tilde )[mss.ensemble.index]




  
  ## derivative of xi.ll w.r.t. f0
  dxi.ll.daf0 = ( dxi.ll.dC * C / af0 ) [mss.ensemble.index]


  dFactor.1.daf0 = f0Koaf0 - fmKoaf0 * xi.ss
  dFactor.1.df0K = af0.i
  dFactor.1.dfmK = af0.i * xi.ss

  dFactor.2.daf0 = ( -0.75 * ( log(xi.ll.i) + 1 )  + (b0K + bmK * xi.ss) ) * dxi.ll.daf0 - 2/af0.i * bmK * xi.ss * xi.ll.i
  dFactor.2.db0K = xi.ll.i
  dFactor.2.dbmK = xi.ss * xi.ll.i

  
  dFactor.3.daf0 = (
                  - 0.75 * xi.ll.times.g.tilde[mss.ensemble.index]/xi.ll.i * dxi.ll.daf0 -
                  0.75 * xi.ll.times.g.tilde.deri[mss.ensemble.index] * (
                                                                      C.i/2/sqrt( xi.ll.i) * dxi.ll.daf0 + sqrt(xi.ll.i) * C.i/af0.i
                                                                      )
                  )
  
  F12 = Factor.1 * Factor.2
  F23 = Factor.2 * Factor.3
  F13 = Factor.1 * Factor.3
  
  dfK.daf0.raw = dFactor.1.daf0 * F23 + dFactor.2.daf0 * F13 + dFactor.3.daf0 * F12
  dfK.df0K = dFactor.1.df0K * F23
  dfK.dfmK = dFactor.1.dfmK * F23
  dfK.db0K = dFactor.2.db0K * F13
  dfK.dbmK = dFactor.2.dbmK * F13


  ## bring the derivative w.r.t. a in to the correct column of the gradient
  dfK.daf0 = array(0, dim = c(length(mss.ensemble.index), aargs$num.ls ))

  for( i in 1:length(mss.ensemble.index) )
    dfK.daf0[i,aargs$ls.index[mss.ensemble.index[i]]] = dfK.daf0.raw[i]
  
  
  
  return(
##         dxi.ll.df0
         cbind(rep(0,length(mssrange) ),dfK.daf0 , dfK.df0K, dfK.dfmK, dfK.db0K, dfK.dbmK)
         )

  
}


fpifK.2 <- function(x,par,aargs){
  cache=fpi.2.calc.cache(x,par,aargs)
  return(c(fpi.2(x,par,aargs,cache),fK.2(x,par,aargs,cache)))
}





dfpifK.2 <- function(x,par,aargs){

  
pionrange = 1:aargs$num.ens
  
  cache=dfpi.2.calc.cache(x,par,aargs)

  dfpifK <- array(0, dim = c(length(x),5+aargs$num.ls ) )

  dfpifK[pionrange,1:(1+aargs$num.ls)] = dfpi.2(x[pionrange],par,aargs,cache)
  kaonrange = ( 1:length(aargs$mss.ensemble.index)) + aargs$num.ens
##  print(kaonrange)
  dfpifK[ kaonrange, ]  = dfK.2(x,par,aargs,cache)
  
  return(dfpifK)
}
