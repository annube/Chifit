library(elliptic,pos = "package:base")

hbarc <- 197.6 ## MeV x fm


f.pi.exp <- 130.7
m.pi.exp <- 135.0

## the following 3 functions do basically the same as the 3 preceding ones but they operate on xi = m_pi^2/ ( 4 pi f_0)^2
## with C = 4 * pi * f_0 * a
xi.FV_fn <- function( xi.IV,C){
  res <-  xi.IV*(1+xi.IV*g.tilde( C * sqrt(xi.IV) ,partitions) )
  return(res)
}

dxi.FV_fn <- function( xi.IV,C){
  xi.FV <-  xi.IV*(1+xi.IV*g.tilde( C * sqrt(xi.IV) ,partitions) )
  return(
         2*xi.FV/xi.IV - 1 +
         C * xi.IV^(1.5) / 2 * g.tilde.deri( C * sqrt(xi.IV),partitions )
         )
}

determine.xi.IV <- function(xi.FV,C){


##    print("----")
##      print(xi.FV)
##      print(C)
##    print("----")

  ## make a good first guess

  xI <- xi.FV

   for( i in 1:2) {
     G = g.tilde( sqrt(xI) * C,partitions)
     xI =   1/ (2 * G) * ( sqrt( 1 + 4 * G * xi.FV ) - 1 )
   }

  ## and then use newton method to find the zero

  res <- numeric(length(xi.FV))

  for( i in 1:length(res) )
    res[i] <- newton_raphson(xI[i],function(x) xi.FV_fn(x,C)-xi.FV[i], function(x) dxi.FV_fn(x,C),100,tol=1.e-10)$root
  
  return(res)
  
}

## derivative of determine.xi.IV w.r.t. C
ddetermine.xi.IV.dC <- function(xi.FV,xi.IV,C,cache=NULL){

  if( !is.null(cache) ) {
    g.tilde.precomp = cache$xi.ll.times.g.tilde / xi.IV
    g.tilde.deri.precomp = cache$xi.ll.times.g.tilde.deri / xi.IV
  } else {
    g.tilde.precomp = ( xi.FV/xi.IV - 1 ) / xi.IV
    print(g.tilde.precomp)
    print( g.tilde( sqrt(xi.IV) * C ,partitions ) )
    g.tilde.deri.precomp = g.tilde.deri( sqrt(xi.IV) * C ,partitions )
  }
    
  return(
         - (
            2/C*xi.FV
            + 
            (xi.IV^(2.5)) * g.tilde.deri.precomp
            ) / (
                 1 +
                 2 * xi.IV*g.tilde.precomp
                 +
                 C/2 * ( xi.IV^(1.5) ) * g.tilde.deri.precomp
                 )
         )
  
}


## chi-PT formula for the pion decay constant
fpi <- function(x,par,aargs,cache=NULL){

  
  nens <- aargs$num.ens
  pionrange <- 1:nens
  
  f0 <- par[1]
  a <- par[aargs$ls.index[pionrange]+1]
  C= 4*pi*f0 *aargs$Loa*a
  xi.ll.FV <- (x[pionrange]/(a[pionrange]*4*pi*f0))^2
  

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

  ## constraint on b
  xi.ll.exp <- (m.pi.exp)^2/(4*pi*f0)^2

  
  b <- ( f.pi.exp / f0 - 1  ) / xi.ll.exp + 2 * log( xi.ll.exp )

                    
  return(
          a*f0 *
          ( 1-2*xi.ll*log(xi.ll) + xi.ll  *b ) *
          (1-2 * xi.ll.times.g.tilde )
##         xi.ll
         )
}

## derivative of fpi w.r.t. parameters par
dfpi <- function(x,par,aargs,cache=NULL){

  
  nens <- aargs$num.ens
  pionrange <- 1:nens
  
  f0 <- par[1]
  a <- par[aargs$ls.index[pionrange]+1]
  C= 4*pi*f0 *aargs$Loa*a[pionrange]
  xi.ll.FV <- (x[pionrange]/(a[pionrange]*4*pi*f0))^2
  
  ## receycle precomputed xi.ll
  if( is.null(cache) ) {
    
    xi.ll <- numeric(length(pionrange))
    for( i in 1:nens ) {
      xi.ll[i] = determine.xi.IV( xi.ll.FV[i] , C[i])
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

  

  
  ## constraint on b
  xi.ll.exp <- (m.pi.exp)^2/(4*pi*f0)^2

  ## constraint from physical values imposed on b
  b <- ( f.pi.exp / f0 - 1  ) / xi.ll.exp + 2 * log( xi.ll.exp )

  ## the 3 factors of which f_PS is being constructed
  Factor.1 = a*f0
  Factor.2 = (1-2*xi.ll*log(xi.ll) + b*xi.ll )
  Factor.3 = (1-2 * xi.ll.times.g.tilde )

  ## derivative of xi.ll w.r.t. f0
  dxi.ll.df0 = dxi.ll.dC * C / f0

  ## derivative of the three factors w.r.t. f_0
  dFactor.1.df0 = a
  dFactor.2.df0 = ( - 2 * ( log( xi.ll) + 1 ) + b ) * dxi.ll.df0 +
    xi.ll * ( 1/f0 * ( 1/xi.ll.exp * ( f.pi.exp/f0-2) - 4 ) )
  dFactor.3.df0 = -2 * xi.ll.times.g.tilde/xi.ll * dxi.ll.df0 -
    2 * xi.ll.times.g.tilde.deri * (
                                    C/2/sqrt( xi.ll) * dxi.ll.df0 + sqrt(xi.ll) * C/f0
                                    )


  ## derivative of xi.ll w.r.t. a
  dxi.ll.da = dxi.ll.dC * C / a

  ## derivative of the three factors w.r.t. a
  dFactor.1.da = f0
  dFactor.2.da = ( - 2 * ( log( xi.ll) + 1 ) + b ) * dxi.ll.da
  dFactor.3.da = -2 * xi.ll.times.g.tilde/xi.ll * dxi.ll.da -
    2 * xi.ll.times.g.tilde.deri * (
                                    C/2/sqrt( xi.ll) * dxi.ll.da + sqrt(xi.ll) * C/a
                                    )

  ## construct complete derivative w.r.t. a
  dfpi.da.raw =     dFactor.1.da  *  Factor.2     *  Factor.3 +
                 Factor.1 * dFactor.2.da *  Factor.3 +
                 Factor.1 *  Factor.2     * dFactor.3.da

  ## bring the derivative w.r.t. a in to the correct column of the gradient
  dfpi.da = array(0, dim = c(nens, aargs$num.ls ))

  for( i in 1:nens )
    dfpi.da[i,aargs$ls.index[i]] = dfpi.da.raw[i]


  ## construct complete derivative w.r.t. f_0  
  dfpi.df0 =    dFactor.1.df0  *  Factor.2     *  Factor.3 +
                 Factor.1 * dFactor.2.df0 *  Factor.3 +
                 Factor.1 *  Factor.2     * dFactor.3.df0


  ## bring the derivative w.r.t. a in to the correct column of the gradient
  dxi.ll.da.s = array(0, dim = c(nens, aargs$num.ls ))

  for( i in 1:nens )
    dxi.ll.da.s[i,aargs$ls.index[i]] = dxi.ll.da[i]

  
                    
  return(
         ## combine to columns
          cbind(dfpi.df0,dfpi.da)
##         cbind(dxi.ll.df0,dxi.ll.da.s)
         )
}


fK <- function(x,par,aargs,cache=NULL){

  nens <- aargs$num.ens
  pionrange <- 1:nens
  mssrange=(1:aargs$num.mss)+nens
  
  f0 <- par[1]
  a <- par[aargs$ls.index[pionrange]+1]

  f0K = par[(1+aargs$num.ls)+1]
  fmK = par[(1+aargs$num.ls)+2]

  b0K = par[(1+aargs$num.ls)+3]
  bmK = par[(1+aargs$num.ls)+4]



  C= 4*pi*f0 *aargs$Loa*a
  xi.ll.FV <- (x[pionrange]/(a*4*pi*f0))^2

  
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

  
  xi.ss <- (x[mssrange]/a[mss.ensemble.index])^2/(4*pi*f0)^2


  Factor.1 = a[mss.ensemble.index] * ( f0K + fmK * xi.ss )
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

dfK <- function(x,par,aargs,cache=NULL){

  nens <- aargs$num.ens
  pionrange <- 1:nens
  mssrange=(1:aargs$num.mss)+nens
  
  f0 <- par[1]
  a <- par[aargs$ls.index[pionrange]+1]

  f0K = par[(1+aargs$num.ls)+1]
  fmK = par[(1+aargs$num.ls)+2]

  b0K = par[(1+aargs$num.ls)+3]
  bmK = par[(1+aargs$num.ls)+4]



  xi.ll.FV <- (x[pionrange]/(a*4*pi*f0))^2
  C= 4*pi*f0 *aargs$Loa*a[pionrange]
  C.i = C[mss.ensemble.index]
  
  if( is.null(cache) ) {
    
    xi.ll <- numeric(length(pionrange))
    for( i in 1:aargs$num.ens ) {
      xi.ll[i] = determine.xi.IV( (x[i]/(a[i]*4*pi*f0))^2 , C[i])
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
  a.i <- a[mss.ensemble.index]
  xi.ss <- (x[mssrange]/a[mss.ensemble.index])^2/(4*pi*f0)^2


  Factor.1 = a[mss.ensemble.index] * ( f0K + fmK * xi.ss )
  Factor.2 = (
              1
              - 0.75*(xi.ll.i*log(xi.ll.i))
              + (b0K+bmK*xi.ss)*xi.ll.i
              )
  Factor.3 = (1-3/4* xi.ll.times.g.tilde )[mss.ensemble.index]




  
  ## derivative of xi.ll w.r.t. f0
  dxi.ll.df0 = ( dxi.ll.dC * C / f0 ) [mss.ensemble.index]
  ## derivative of xi.ll w.r.t. a
  dxi.ll.da = ( dxi.ll.dC * C / a) [mss.ensemble.index]


  dFactor.1.df0 =     - a.i * 2 * fmK / f0* xi.ss
  dFactor.1.da = f0K - fmK * xi.ss
  dFactor.1.df0K = a.i
  dFactor.1.dfmK = a.i * xi.ss

  dFactor.2.df0 = -0.75 * ( log(xi.ll.i) + 1 ) * dxi.ll.df0 - 2/f0 * bmK * xi.ss * xi.ll.i + (b0K + bmK * xi.ss) * dxi.ll.df0
  dFactor.2.da = -0.75 * ( log(xi.ll.i) + 1 ) * dxi.ll.da - 2/a.i * bmK * xi.ss * xi.ll.i + (b0K + bmK * xi.ss) * dxi.ll.da
  dFactor.2.db0K = xi.ll.i
  dFactor.2.dbmK = xi.ss * xi.ll.i

  
  dFactor.3.da = (
                  - 0.75 * xi.ll.times.g.tilde[mss.ensemble.index]/xi.ll.i * dxi.ll.da -
                  0.75 * xi.ll.times.g.tilde.deri[mss.ensemble.index] * (
                                                                      C.i/2/sqrt( xi.ll.i) * dxi.ll.da + sqrt(xi.ll.i) * C.i/a.i
                                                                      )
                  )
  dFactor.3.df0 = (
                  - 0.75 * xi.ll.times.g.tilde[mss.ensemble.index]/xi.ll.i * dxi.ll.df0 -
                  0.75 * xi.ll.times.g.tilde.deri[mss.ensemble.index] * (
                                                                      C.i/2/sqrt( xi.ll.i) * dxi.ll.df0 + sqrt(xi.ll.i) * C.i/f0
                                                                      )
                  )
  F12 = Factor.1 * Factor.2
  F23 = Factor.2 * Factor.3
  F13 = Factor.1 * Factor.3
  
  dfK.df0 = dFactor.1.df0 * F23 + dFactor.2.df0 * F13 + dFactor.3.df0 * F12
  dfK.da.raw = dFactor.1.da * F23 + dFactor.2.da * F13 + dFactor.3.da * F12
  dfK.df0K = dFactor.1.df0K * F23
  dfK.dfmK = dFactor.1.dfmK * F23
  dfK.db0K = dFactor.2.db0K * F13
  dfK.dbmK = dFactor.2.dbmK * F13


  ## bring the derivative w.r.t. a in to the correct column of the gradient
  dfK.da = array(0, dim = c(length(mss.ensemble.index), aargs$num.ls ))

  for( i in 1:length(mss.ensemble.index) )
    dfK.da[i,aargs$ls.index[mss.ensemble.index[i]]] = dfK.da.raw[i]
  
  
  
  return(
##         dxi.ll.df0
         cbind(dfK.df0 , dfK.da, dfK.df0K, dfK.dfmK, dfK.db0K, dfK.dbmK)
         )

  
}

dfpifK <- function(x,par,aargs){

  
  pionrange <- 1:aargs$num.ens

  f0 <- par[1]
  a <- par[aargs$ls.index[pionrange]+1]


  C= 4*pi*f0 *aargs$Loa*a[pionrange]

  xi.ll.FV = (x[pionrange]/(a[pionrange]*4*pi*f0))^2
  
  xi.ll <- numeric(length(pionrange))
  for( i in 1:aargs$num.ens ) {
    xi.ll[i] = determine.xi.IV( xi.ll.FV[i] , C[i])
  }
  
  xi.ll.times.g.tilde = ( xi.ll.FV/xi.ll - 1 )  ## xi.ll * g.tilde(sqrt(xi.ll)*C,partitions)
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

  dfpifK <- array(0, dim = c(length(x),5+aargs$num.ls ) )

  dfpifK[pionrange,1:(1+aargs$num.ls)] = dfpi(x[pionrange],par,aargs,cache)
  kaonrange = ( 1:length(aargs$mss.ensemble.index)) + aargs$num.ens
##  print(kaonrange)
  dfpifK[ kaonrange, ]  = dfK(x,par,aargs,cache)
  
  return(dfpifK)
}


fpifK <- function(x,par,aargs){

  
  pionrange <- 1:aargs$num.ens

  f0 <- par[1]
  a <- par[aargs$ls.index[pionrange]+1]


  C= 4*pi*f0 *aargs$Loa*a

  xi.ll <- numeric(length(pionrange))
  for( i in 1:aargs$num.ens ) {
    xi.ll[i] = determine.xi.IV( (x[i]/(a[i]*4*pi*f0))^2 , C[i])
  }
  
  xi.ll.times.g.tilde = xi.ll*g.tilde(sqrt(xi.ll)*C,partitions)
  
  
  cache=list( xi.ll = xi.ll ,
    xi.ll.times.g.tilde = xi.ll.times.g.tilde 
    )

  
  return(c(fpi(x,par,aargs,cache),fK(x,par,aargs,cache)))
}



dfK_ <- function(x,par,aargs){
  dfK(c(aargs$mpi,x),par,aargs)
}

fK_ <- function(x,par,aargs){
  fK(c(aargs$mpi,x),par,aargs)
}
