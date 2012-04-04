
cosh.R <- function(par,t,T) cosh( par[1] * (t - T/2 ) ) /  cosh( par[1] * (t + 1 - T/2 ) )
dcosh.R <- function(par,t,T) {
  (
   (t - T/2) * sinh( par[1] * (t - T/2 ) ) * cosh( par[1] * (t + 1 - T/2 ) ) -
   
   cosh( par[1] * (t - T/2 ) ) * (t + 1 - T/2 ) *  sinh( par[1] * (t + 1 - T/2 ) )
   
   )  /  cosh( par[1] * (t + 1 - T/2 ) )^2
}


fit.cosh.approx <- function(corr,dcorr,range,T){
  res = numeric(length(range)-1)
  for( i in 1:(length(range)-1) ){
    t = (0:(T/2))[range[i]]
    R = corr[range[i]]/corr[range[i+1]]
    res[i] <- newton(function(x,...) cosh.R(x,...) - R , function(x,...) dcosh.R(x,...) ,1,t=t,T=T)
  }
  return(res);
}


corr.cosh <- function(x,par,T){
  par[1] * ( exp( - par[2] * x) + exp( - par[2] * (T-x) ) )
}


dcorr.cosh <- function(x,par,T) {
cbind(
                         exp( - par[2] * x) + exp( - par[2] * (T-x) ) ,
                         par[1] * ( -x * exp( - par[2] * x) + -( T - x ) * exp( - par[2] * (T-x) ) )
                         )
}

corr.exp <- function(x,par,...) par[1] * exp(-par[2]*x)


dcorr.exp <- function(x,par,...) cbind( exp(-par[2]*x), - x* par[1] * exp(-par[2]*x) ) 


make.corr.fit <- function(corr,dcorr,T=2*(length(corr)-1),from.T=T/4,to.T=0.4*T,fntype="cosh"){
  range=min(from.T,to.T):max(from.T,to.T)

  xs = 0:(T/2)

  len.ra = length(range)
  
  if( fntype == "cosh" ){
    f = function(x,par,...) corr.cosh(x,par,T=T)
    df = function(x,par,...) dcorr.cosh(x,par,T=T)
    meffs <- fit.cosh.approx(corr,dcorr,range,T)
    meff.guess = mean(meffs)
    par1.guess = mean( (corr[range])/ corr.cosh(xs[range],c(1,meff.guess),T=T) )
  } else {
    f = corr.exp
    df = dcorr.exp
    meffs <- log(corr[range[1:(len.ra-1)]]/corr[range[2:(len.ra)]])
    meff.guess = mean(meffs)
    par1.guess = mean( (corr[range])/ corr.exp(xs[range],c(1,meff.guess)) )
    ##par1.guess = 1
  }

##  print( meff.guess)
##  print( par1.guess)

  if(length(range) > 2 ) {
    wnlls.res <- wnlls(
                       xs
                       ,
                       corr,
                       dy=dcorr,
                       f= f,
                       df= df ,
                       par=c(par1.guess,meff.guess),
                       range=range
                       )
    return( wnlls.res )
  } else {
    return(NA)
  }
  
}


determine.optimal.fit.range <- function(corr,dcorr,T,min.conf = 0.1,tmin=2,p = 0.95){


  tmax = max( which( abs( dcorr / corr ) < min.conf ) )
##  print(tmax)
  
  res.all.chisqr <- array(NA,dim=c(T/2+1,T/2+1))
  chisqr.test = array(NA,dim=c(T/2+1,T/2+1))
  res.all.compat.zero = array(NA,dim=c(T/2+1,T/2+1))
  valid.i.index = array(NA,dim=c(T/2+1,T/2+1))
  valid.j.index = array(NA, dim = c(T/2+1,T/2+1)) 

  valid.i.index.just.chi = array(NA,dim=c(T/2+1,T/2+1))
  valid.j.index.just.chi = array(NA, dim = c(T/2+1,T/2+1)) 

  for( i in tmin:tmax) {
    if( tmax >= i+3)
      for(j in (i+3):tmax ){
        res <- try( make.corr.fit(corr,dcorr,T,i,j) )
        if( ! is.na(res) &  ! inherits( res,"try-error") ) {
          res.all.chisqr[i,j] = res$Chisqr
          chisqr.test[i,j] =  pchisq( res$Chisqr , df = j-i+1-2 ) < p
          range = i:j
          res.all.compat.zero[i,j] = all( abs( (corr-res$predict)[range]/dcorr[range] ) < 1.96 )
##          print( paste( i , " " , j , " " , res$Chisqr ) )
          if( res.all.compat.zero[i,j] & chisqr.test[i,j] ) {
            valid.i.index [i,j] = i
            valid.j.index [i,j] = j
          }
          if( chisqr.test[i,j] ) {
##             print(paste("i " , i," j ", j))
            valid.i.index.just.chi [i,j] = i
            valid.j.index.just.chi [i,j] = j
          }
        }
      }
  }
  non.na.range = ! is.na(valid.i.index)
  ij.combinations <-  cbind( valid.i.index[non.na.range] , valid.j.index[non.na.range] )

  
  non.na.range = ! is.na(valid.i.index.just.chi)
  ij.combinations.just.chi <-  cbind( valid.i.index.just.chi[non.na.range] , valid.j.index.just.chi[non.na.range] )

  ij.diff = ij.combinations[,2]-ij.combinations[,1]
  max.dof.index <- which( max(ij.diff) == ij.diff )

  return(list( ij.combinations = ij.combinations,
              ij.combinations.just.chi = ij.combinations.just.chi,
              max.dof.index = max.dof.index
              )
         )

}
