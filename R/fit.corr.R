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


  n.sim = floor( dim(data)[1] / boot.l ) * boot.l
    
  ts.boot.res <- tsboot(data,
                        boot.mean ,
                        boot.R ,
                        sim = "fixed" ,
                        l = boot.l ,
                        parallel = "multicore",
                        ncpus=ncpus,
                        n.sim = n.sim ,
                        endcorr=FALSE )

  boot.err <- apply( ts.boot.res$t,2,sd)

  
  return( list( std.err = std.err ,
               boot.err = boot.err,
               correction.factor = sqrt( n.sim ) / sqrt( dim(data)[1] ),
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


corr.make.range <- function(t1,t2)
  (t1+1):(t2+1)


corr.mk.range.all <- function( range , ncorr , T )
  rep( range , ncorr ) + rep( 0:(ncorr-1) * (T/2+1) , each = length(range) )




simple.correlator.fit <- function(  Corr , dCorr  , f ,df , T , p.min = 0.1){
  Thalf = T/2

  range = corr.make.range(  round(Thalf/2) ,  Thalf )
  range.all = corr.mk.range.all( range, length( Corr ) , T )
  
  
  lm.res.q <- lev.marq(
                       ( 0:Thalf ) [ range ] ,
                       unlist( Corr )[range.all] ,
                       unlist ( dCorr ) [range.all] ,
                       f ,df,
                       c( rep( 0.1 , length(Corr) ) , 0.1)
                       )


  p.value <- 0
  
  t1 = 0

  while( p.value < p.min ) {
    t1 = t1 + 1
    range = corr.make.range(  t1 ,  Thalf )
    range.all = corr.mk.range.all( range, length( Corr ) , T )
    lm.res <- lev.marq( ( 0:Thalf ) [ range ] , unlist( Corr )[range.all] , unlist ( dCorr ) [range.all] , f ,df,
                       lm.res.q$beta )

    p.value = 1 - pchisq( lm.res$Chisqr , length( range.all ) - length( lm.res$beta ) )
##    print( p.value )
    
  }
  
  
  return( list( lm.res = lm.res , t1 = t1 ) )
  
}



fit.correlator <- function( data  , T ,t1 , t2 , t1.exc1 = round(t1/2), automatic.t1.exc1.adjust = TRUE , max.chisqr.o.dof = 0.05, ncpus = 4 , max.exc.state.contamin = 0.2 ,...) {
  require(hadron)

  if( length(data) == 1  ){
    mycosh.fn = mycosh
    dmycosh.fn = dmycosh
    num.par = 2
  } else if( length(data) == 2 ) {
    mycosh.fn = mycosh.2
    dmycosh.fn = dmycosh.2
    num.par = 3
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
######   step  : construct two states correlation function and fit parameters
######             of ground state and first excited state simultaneously
######
######
  
  
  mycosh.two.states <- function( x, par )
    mycosh.fn(x,par[1:num.par],list(T=T)) + mycosh.fn(x,par[(1:num.par) + num.par],list(T=T))

  dmycosh.two.states <- function( x, par )
    cbind(dmycosh.fn(x,par[1:num.par],list(T=T)) , dmycosh.fn(x,par[(1:num.par)+num.par],list(T=T)) )



  corr.all <- unlist ( corr )
  dcorr.all <- unlist ( dcorr )

  scf.res <- simple.correlator.fit( corr,dcorr ,
                                   function(x,par)  mycosh.fn( x, par , list(T=T) ) , 
                                   function(x,par)  dmycosh.fn( x, par , list(T=T) ) ,
                                   T , p.min = 0.95 )
  
  
  fit.lm.two.states <- function( corr.data , t1.fit = t1.exc1 , t2.fit = T/2 , 
                                par0 = c(  scf.res$lm.res$beta ,  scf.res$lm.res$beta+0.5 ) )
    {
      
      range = corr.make.range( t1.fit , t2.fit )
      range.all = corr.mk.range.all( range , length(corr) , T )
      
      lm.res <- lev.marq(  ( 0:(T/2) ) [ range ]  , corr.data[range.all] ,
                         dcorr.all [ range.all ] ,
                         function(x,par) mycosh.two.states(x,par),
                         function(x,par) dmycosh.two.states(x,par),
                         par0 )
    
      return( lm.res )
  }



  ## perform one fit with the initial value of t1.exc1 passed to the function to determine
  ## initial set of parameters

  lm.res.two.states.i <- fit.lm.two.states(corr.all,t1.exc1 , T/2)


  if ( automatic.t1.exc1.adjust ) {
  
    p.value = 0
    chisqr.o.dof = 1
    t1.exc1 = 0
    
    while( chisqr.o.dof > max.chisqr.o.dof ){

      t1.exc1 = t1.exc1 + 1
      lm.res.two.states <- fit.lm.two.states(corr.all,t1.exc1 , T/2 , lm.res.two.states.i$beta )

      dof = length(corr) * ( T/2 - t1.exc1 + 1 ) - length( lm.res.two.states$beta )
      p.value =
        1  -  pchisq(
                     lm.res.two.states$Chisqr ,
                     dof
                     )

      chisqr.o.dof = lm.res.two.states$Chisqr / dof
      
    }
    print( paste( "found optimal t1.exc1 = " , t1.exc1 ) )
  }
  else {
    lm.res.two.states <- fit.lm.two.states(corr.all,t1.exc1 , T/2 , lm.res.two.states.i$beta )
  }

  
  dof.two.states = (T/2-t1.exc1 + 1 ) * length(data) - length( lm.res.two.states$beta )
  

  ## prepare table of boot strap sampled correlators
  
  corr.all.boot <- array( 0 , dim=c( boot.R , (T/2 + 1 ) * length(data) ) )


  for( lei in 1:length(data) ) {
    corr.all.boot[, (lei-1) * (T/2+1) + 1:(T/2+1)] = ce.res[[lei]]$boot.res$t
  }

  
  lm.res.two.states.boot <- mcapply(
    corr.all.boot ,
    function( corr ) { fr <- fit.lm.two.states(corr , par0 =  lm.res.two.states$beta); return( c( fr$beta,fr$Chisqr ) ) } , mc.cores=ncpus, mc.preschedule = FALSE
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

  range = corr.make.range( t1.exc1 , T/2 )

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




  


######
######
######
######   step 2a : find optimal t1 for definition of a region without contamination of the excited state
######
######


  par.range.exc1 = (1+length(data)) + 1:( (1+length(data)) )


  t1.opt  = numeric( length(data) )

  for( lei in 1:length(data) ) {
    predict.exc1 = mycosh( 0:(T/2) , ( lm.res.two.states$beta[par.range.exc1])[c(lei,length(data)+1)] , list(T=T) )

    t1.opt[lei] = min( which(   ( predict.exc1 / dcorr[[lei]] < max.exc.state.contamin ) & ( 0:(T/2) > t1.exc1 ) ) )  - 1
  }

  if( any( is.infinite(t1.opt) ) )
    t1.opt[] = T/2-5
  
  print( t1.opt )


  t1 = max( t1.opt )

    
  
######
######
######
######   step 2 : fit correlator directly in ground state region ##########
######
######


  

  range = corr.make.range ( t1 , t2  )



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

  corr.all.range.boot <- array( 0 , dim=c( boot.R , length(range) * length(data) ) )


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
  

  
  return(
         list(
              ce.res = ce.res,
              lm.res = lm.res,
              lm.res.2 = lm.res.2,
              lm.res.boot = lm.res.boot,
              dof = dof,
              t1 = t1 , t1.exc1 = t1.exc1, t2 = t2,
              lm.res.two.states = lm.res.two.states,
              lm.res.two.states.boot = lm.res.two.states.boot,
              dof.two.states = dof.two.states
              )
         )
  
  
}

