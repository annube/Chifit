## the final fuction for fitting correlators


mycosh <- function(x,par,aargs) {
  0.5 * par[1] * ( exp( - par[2] * x ) + exp( - par[2] * (aargs$T-x) ) )
}

dmycosh <- function(x,par,aargs) {
  cbind(    0.5* ( exp( - par[2] * x ) + exp( - par[2] * (aargs$T-x) ) ) ,
        0.5 * par[1] * ( - x *  exp( - par[2] * x ) - (aargs$T -  x) * exp( - par[2] * (aargs$T-x) ) )
        )
}



correlator.error <- function( data , boot.R = 100 , boot.l = 10 ) {
  require(boot)

  std.err <- apply( data , 2 , sd ) / sqrt( dim(data)[1] )
  
  boot.mean <- function(data)
    apply( data , 2 , mean )

    
  ts.boot.res <- tsboot(data, boot.mean , boot.R , sim = "fixed" , l = boot.l ,parallel = "multicore",ncpus=4 )

  boot.err <- apply( ts.boot.res$t,2,sd)

  
  return( list( std.err = std.err ,
               boot.err = boot.err,
               boot.res = ts.boot.res
               )
         )  
  
}


mcapply <- function( X , fun , ... ) {
  res.l <- mclapply( 1:dim(X)[1] , function( i ) fun( X[i,] ), ... )
  
  res <- array( dim = c(length(res.l),length(res.l[[1]])) )
  for( i in 1:dim(res)[1] ) res[i,] = res.l[[i]]
  return( res )
}


fit.correlator <- function( data , T ,t1 , t2 , num.t.points = 4,correlated.fit = T, ...) {
  require(multicore)
  require(hadron)

  corr <- apply( data ,2 ,mean)
  ce.res <- correlator.error( data , ...)
  dcorr <- ce.res$boot.err

  range.mo = (t1+1):(t2)
  
  m.approx <- fit.cosh.approx( corr,dcorr,1:dim(data)[2], T )

  m.approx.boot <- t(  mcapply( ce.res$boot.res$t , function( corr )
                            fit.cosh.approx( corr ,dcorr,1:dim(data)[2], T ) ,
                            mc.cores = 4 , mc.preschedule = FALSE 
                            )
                     )

  dm.approx <- apply( m.approx.boot,1,sd)

  m.approx.mean <- lm(m.approx[range.mo]~1,weights=1/dm.approx[range.mo]^2)$coefficients[1]

  chisqr = sum( ( ( m.approx - m.approx.mean ) / dm.approx )[range.mo]^2 )
  
  if( 1-pchisq( chisqr , length(range.mo) - 1 ) < 0.05 ) {
    print( "Waring p-value of meff fit is < 0.05 . This usually means that the T-range should be adjusted " )
    print( chisqr )
    print( length(range.mo)-1)
  }

  plotwitherror( range.mo-0.5 , m.approx[range.mo] , dm.approx[range.mo] )
  abline( h = m.approx.mean )


  range = round( seq( (t1+1),(t2+1),length.out = num.t.points ) )

  
  lm.res <- lev.marq( ( 0:(T/2+1) )[range] , corr[range],dcorr[range],
                     function(x,par) mycosh(x,par,list( T = T )),
                     function(x,par) dmycosh(x,par,list( T = T )),
                     c( 1, 1 )
                     )

  
  if( correlated.fit ) {
    C <- cov( ce.res$boot.res$t )
    print( sqrt(diag(C)) / dcorr )
  }
  else
    C <- diag( dcorr^2 )

  for( i in 1:2 ){
  
    fit.wlm <- function( corr.data ) {
      wlm.res <- wnlls(  0:(T/2)  , corr.data,
                       C = C,
                       f.in = function(x,par) mycosh(x,par,list( T = T )),
                       df.in = function(x,par) dmycosh(x,par,list( T = T )),
                       par = lm.res$beta,
                       range = range
                       )
      return( wlm.res )
    }
    
    wlm.res <- fit.wlm( corr )
   
    wlm.res.boot <- mcapply( ce.res$boot.res$t , function( corr ) { fr <- fit.wlm(corr); return( c( fr$beta,fr$Chisqr ) ) } , mc.cores=4, mc.preschedule = FALSE)
    
  }
  

  plot.wnlls(wlm.res,log='y')
  
  return(
    list(
      ce.res = ce.res,
      m.approx = m.approx, m.approx.boot = m.approx.boot ,
      wlm.res = wlm.res,
      wlm.res.boot = wlm.res.boot
      )
    )
  
  
}

