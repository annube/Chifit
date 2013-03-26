## the final fuction for fitting correlators


mycosh <- function(x,par,aargs) {
  0.5 * par[1] * ( exp( - par[2] * x ) + exp( - par[2] * (aargs$T-x) ) )
}

dmycosh <- function(x,par,aargs) {
  cbind(    0.5* ( exp( - par[2] * x ) + exp( - par[2] * (aargs$T-x) ) ) ,
        0.5 * par[1] * ( - x *  exp( - par[2] * x ) - (aargs$T -  x) * exp( - par[2] * (aargs$T-x) ) )
        )
}

mycosh.2 <- function(x,par,aargs) {
 c( mycosh( x, c( par[1] , par[3] ) , aargs ) , mycosh( x, c( par[2] , par[3] ) , aargs ) )
}

dmycosh.2 <- function(x,par,aargs) {
  res <- matrix( 0 , nrow = 2 * length(x) , ncol = 3 )
  one.xrange = 1:length(x)
  df1 = dmycosh( x, c( par[1] , par[3] ) , aargs )
  df2 = dmycosh( x, c( par[2] , par[3] ) , aargs )
  res[ one.xrange , 1] = df1[,1]
  res[ one.xrange , 3] = df1[,2]
  res[ one.xrange + length(x) , 2] = df2[,1]
  res[ one.xrange + length(x) , 3] = df2[,2]
  return( res )
}



correlator.error <- function( data , boot.R = 100 , boot.l = 10,ncpus=4 ) {
  require(boot)

  std.err <- apply( data , 2 , sd ) / sqrt( dim(data)[1] )
  
  boot.mean <- function(data)
    apply( data , 2 , mean )

    
  ts.boot.res <- tsboot(data, boot.mean , boot.R , sim = "fixed" , l = boot.l ,parallel = "multicore",ncpus=ncpus )

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


fit.correlator <- function( data  , T ,t1 , t2 , t1.exc1 = round(t1/2), automatic.t1.adjust = rep(TRUE,2) , num.t.points = 4,correlated.fit = T, ncpus = 4 ,...) {
  require(hadron)
  
  if( length(data) == 1  ){
    mycosh.fn = mycosh
    dmycosh.fn = dmycosh
  } else if( length(data) == 2 ) {
    mycosh.fn = mycosh.2
    dmycosh.fn = dmycosh.2
  }
    
  
######
######
######
######   step 0 : calculate correlator error ##########
######
######
  
  
  corr <- lapply( data ,function(le) apply( le , 2 ,mean) )
  ce.res <- lapply( data , function( le ) correlator.error( le , ncpus = ncpus, ...) )
  boot.R <-  ce.res[[1]]$boot.res$R
  dcorr <- lapply( ce.res , function( le ) le$boot.err )



######
######
######
######   step 1 : calculate/fit effective mass of ground state directly and find a plateau ##########
######
######
  
  m.approx = list()
  m.approx.boot = list()
  dm.approx = list()
  
  for( i in 1:length(data) ) {
  
    m.approx[[i]] <- fit.cosh.approx( corr[[i]],dcorr[[i]],1:dim(data[[i]])[2], T )

    
    m.approx.boot[[i]] <- t(  mcapply( ce.res[[i]]$boot.res$t , function( corr )
                                 fit.cosh.approx( corr ,dcorr[[i]],1:dim(data[[i]])[2], T ) ,
                                 mc.cores = ncpus , mc.preschedule = FALSE 
                                 )
                       )
    
    dm.approx[[i]] <- apply( m.approx.boot[[i]],1,sd)
  }


  if( automatic.t1.adjust[1] ) {
    t1 = 0
  }

  p.value = 0
  once = TRUE
  while( (p.value<0.05 && automatic.t1.adjust[1] ) || once  ) {
    once = FALSE
    if ( automatic.t1.adjust[1] ) t1 = t1 + 1
    range.mo = (t1+1):(t2)
    
    m.approx.all.range = unlist( lapply( m.approx , function(le) le[range.mo] ) )
    dm.approx.all.range = unlist( lapply( dm.approx , function(le) le[range.mo] ) )
    
    m.approx.mean <- lm(  m.approx.all.range ~1,weights=1/dm.approx.all.range^2)$coefficients[1]
    chisqr = sum( ( (m.approx.all.range - m.approx.mean )/dm.approx.all.range )^2 )
    p.value = 1 - pchisq( chisqr , length(m.approx.all.range) - 1 )
  }




  ## make a plot
  plotwitherror( rep( 1 : ( T / 2 ) - 0.5 , length(data) ) , unlist( m.approx ) , unlist( dm.approx  ),
                main = "Direct fit of m_eff of lowest state" ,
                xlab = expression(x[0]/a) , ylab = expression( m[eff] / a )
                )

  
  plotwitherror( rep( range.mo - 0.5 , length(data) ) , m.approx.all.range , dm.approx.all.range , col = "orange" , rep = TRUE)
  abline( h = m.approx.mean )

  abline ( v = t1 , lty = 2 ,col ="gray" )


######
######
######
######   step 2 : fit correlator directly in ground state region ##########
######
######

  

##  range = round( seq( (t1+1),(t2+1),length.out = num.t.points ) )
  range = ( t1 + 1 ) :  ( t2 + 1 )



  corr.all.range <- unlist ( lapply( corr , function( le ) le[range] ) )
  dcorr.all.range <- unlist ( lapply( dcorr , function( le ) le[range] ) )


  
  lm.res <- lev.marq(  ( 0:(T/2+1) )[range] , corr.all.range,dcorr.all.range,
                     function(x,par) mycosh.fn(x,par,list( T = T )),
                     function(x,par) dmycosh.fn(x,par,list( T = T )),
                     c( rep(1,length(data) ), 1 )
                     )



  corr.all <- unlist ( lapply( corr , function( le ) le ) )
  dcorr.all <- unlist ( lapply( dcorr , function( le ) le ) )
  
  plot.new()
  plot.window( xlim = c( 0 , T/2 ) , ylim = c( min( corr.all ), max( corr.all ) ) , log = "y"  )
  axis(1)
  axis(2)

  title(main =  "Correlator fit "  , xlab = expression( x[0]/a) , ylab = expression( a^6 * C( x[0]/a ) ) )

  for( i in 1:length(data) ) {
    plotwitherror( 0:(T/2) , corr[[i]] , dcorr[[i]]  ,rep = TRUE )
    plot( function(x) mycosh.fn(x , lm.res$beta,list(T=T) )[length(x)*(i-1)  + 1: length(x)] , xlim=c(0,T/2) , add = TRUE )
    plotwitherror( ( 0:(T/2) ) [range] , (corr[[i]])[range] , (dcorr[[i]])[range] , col = "orange" , rep=TRUE )
  }


  abline(v  = t1-0.5 , lty = 2 , col="gray" )
  


######
######
######
######   step 2a : corresponding boot analyses  ##########
######
######

  
  fit.wlm <- function( corr.data ) {
    lm.res <- lev.marq(  ( 0:(T/2) ) [ range ]  , corr.data ,
                       dcorr.all.range  ,
                       function(x,par) mycosh.fn(x,par,list( T = T )),
                       function(x,par) dmycosh.fn(x,par,list( T = T )),
                       lm.res$beta
                       )
    return( lm.res )
  }
  
  lm.res.2 <- fit.wlm( corr.all.range )

  dof = length( range ) * length(data) - length( lm.res$beta )

##  return( lm.res.2 )

  corr.all.range.boot <- array( NaN , dim=c( boot.R , length(range) * length(data) ) )


  for( lei in 1:length(data) ) {
    corr.all.range.boot[, (lei-1) * length(range) + 1:length(range)] = ce.res[[lei]]$boot.res$t[,range]
  }
  
  lm.res.boot <- mcapply(
    corr.all.range.boot ,
    function( corr ) { fr <- fit.wlm(corr); return( c( fr$beta,fr$Chisqr ) ) } , mc.cores=ncpus, mc.preschedule = FALSE
    )
    

#####
#####
#####  ad 2) plot residues of the correlator fit
#####
#####


  lm.res.predict = mycosh.fn( 0:(T/2) , lm.res$beta, list(T=T) )


  plot( rep( 0:(T/2) , length(data) ) , (corr.all - lm.res.predict)/dcorr.all , pch = 3,ylim = c(-10,10) )
  abline( h = c(-1,1 ) )
  abline( h = c(-2,2 ) ,lty=2 )
  abline(v  = t1-0.5 , lty = 2 , col="gray" )
  

######
######
######
######   step 3a : subtract ground state contribution to correlator
######           and fit first excited state  
######
######

  corr.subtract <- corr.all - lm.res.predict
  


  if ( automatic.t1.adjust[2] ) {
  
    t1.exc1 = 0
    p.value = 0

  }

  p.value = 0
  once = TRUE
  while( ( p.value < 0.05 && automatic.t1.adjust[2] )  || once ) {
    once = FALSE
    if ( automatic.t1.adjust[2] ) t1.exc1 = t1.exc1 + 1
    ##range = round( seq( (t1.exc1+1),(t2+1),length.out = num.t.points ) )
    range = ( t1.exc1 + 1 ):( t2 + 1 )
    range.all = rep( range , length(data) ) + rep( 0:(length(data)-1) * (T/2+1) , each = length(range) )
    lm.res.exc1 <- lev.marq( (0:(T/2))[range] ,
                            corr.subtract[range.all] ,
                            dcorr.all[range.all],
                            function(x,par) mycosh.fn(x,par,list(T=T) ) ,
                            function(x,par) dmycosh.fn(x,par,list(T=T) ),
                            lm.res$beta + c(rep(0,length(data)),0.1)
                            )
    p.value = 1 - pchisq( lm.res.exc1$Chisqr , length(range) - length(lm.res.exc1$beta ) )
  }



######
######
######
######   step 3b : construct two states correlation function and fit parameters
######             of ground state and first excited state simultaneously
######
######
  
  
  mycosh.two.states <- function( x, par )
    mycosh.fn(x,par[1:3],list(T=T)) + mycosh.fn(x,par[4:6],list(T=T))

  dmycosh.two.states <- function( x, par )
    cbind(dmycosh.fn(x,par[1:3],list(T=T)) , dmycosh.fn(x,par[4:6],list(T=T)) )


  
  fit.lm.two.states <- function( corr.data  ) {
    lm.res <- lev.marq(  ( 0:(T/2) ) [ range ]  , corr.data ,
                       dcorr.all [ range.all ] ,
                       function(x,par) mycosh.two.states(x,par),
                       function(x,par) dmycosh.two.states(x,par),
                       c(lm.res$beta, lm.res.exc1$beta )
                       )
    return( lm.res )
  }

  lm.res.two.states <- fit.lm.two.states(corr.all[range.all])
##  lm.res.exc1 = lm.res.two.states


##   par.only.excited = lm.res.two.states$beta
##   par.only.excited[1:length(data)] = 0

##   corr.diff.exc = corr.all - mycosh.two.states( 0:(T/2) , par.only.excited 
  
  

  dof.two.states = length( range ) * length(data) - length( lm.res.two.states$beta )
  

  ## prepare table of boot strap sampled correlators
  
  corr.all.range.boot <- array( NaN , dim=c( boot.R , length(range) * length(data) ) )


  for( lei in 1:length(data) ) {
    corr.all.range.boot[, (lei-1) * length(range) + 1:length(range)] = ce.res[[lei]]$boot.res$t[,range]
  }

  
  lm.res.two.states.boot <- mcapply(
    corr.all.range.boot ,
    function( corr ) { fr <- fit.lm.two.states(corr); return( c( fr$beta,fr$Chisqr ) ) } , mc.cores=ncpus, mc.preschedule = FALSE
    )

####
####
#### make a plot of the two states correlator fit
####
####

  plot.new()
  plot.window( xlim = c( 0 , T/2 ) , ylim = c( min( corr.all ), max( corr.all ) ) , log = "y"  )
  axis(1)
  axis(2)

  title(main =  "Correlator fit - two states "  , xlab = expression( x[0]/a) , ylab = expression( a^6 * C( x[0]/a ) ) )

  for( i in 1:length(data) ) {
    plotwitherror( 0:(T/2) , corr[[i]] , dcorr[[i]]  ,rep = TRUE )
    plot( function(x) mycosh.two.states(x , lm.res.two.states$beta)[length(x)*(i-1)  + 1: length(x)] , xlim=c(0,T/2) , add = TRUE )
    plotwitherror( ( 0:(T/2) ) [range] , (corr[[i]])[range] , (dcorr[[i]])[range] , col = "orange" , rep=TRUE )
  }
  
  abline(v  = t1.exc1-0.5 , lty = 2 , col="gray" )
  

####
####
#### make a plot of the two residue of the fit
####
####
  
  lm.res.two.states.predict <- mycosh.two.states( 0:(T/2) , lm.res.two.states$beta )

  plot( rep( 0:(T/2) , length(data) ) , (corr.all - lm.res.two.states.predict)/dcorr.all , pch = 3 , ylim = c(-10,10) )
  abline( h = c(-1,1 ) )
  abline( h = c(-2,2 ) ,lty=2 )

  abline(v  = t1.exc1-0.5 , lty = 2 , col="gray" )

  
  return(
         list(
              ce.res = ce.res,
              m.approx = m.approx, m.approx.boot = m.approx.boot ,
              lm.res = lm.res,
              lm.res.2 = lm.res.2,
              lm.res.boot = lm.res.boot,
              dof = dof,
              lm.res.exc1 = lm.res.exc1,
              t1 = t1 , t1.exc1 = t1.exc1, t2 = t2,
              lm.res.two.states = lm.res.two.states,
              lm.res.two.states.boot = lm.res.two.states.boot,
              dof.two.states = dof.two.states
              )
         )
  
  
}

