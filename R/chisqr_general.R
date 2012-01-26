



chisqr_genral <- function (hist,f,doplot=FALSE){


  n.breaks = length(hist$breaks)
  
  theoretical.counts <- numeric(length(hist$counts))

  theoretical.counts = ( sapply( hist$breaks[2:n.breaks] , function (x) f(x)) - 
    sapply( hist$breaks[1:(n.breaks-1)] , function (x) f(x)) ) * sum(hist$counts)


  range = theoretical.counts >= 5

  if(doplot){
    ylim = c(0,max( c( theoretical.counts,hist$counts)) )
    plot(hist,ylim=ylim)
    points(hist$mids[range],theoretical.counts[range])
    points(hist$mids[!range],theoretical.counts[!range],pch=4)
  }

  chisqr = sum( ( ( hist$counts - theoretical.counts ) ^ 2 / theoretical.counts )[range] )

  return(list(chisqr=chisqr,n.intervals=sum(range),theoretical.counts = theoretical.counts ))

}


chisqr_norm <- function(data,doplot=FALSE){
  hist.res <- hist(data,plot=FALSE)
  ch.res <- chisqr_genral(hist.res, function(x) pnorm(x,mean=mean(data),sd=sd(data)),doplot )

  return( list( chisqr= ch.res$chisqr , dof = ch.res$n.intervals - 2 ) )
}


chisqr_pois <- function(data,interval,n.split=100,doplot=FALSE){

  if(min(data) < interval[1] | max(data) > interval[2] ) {
    print(" error in chisqr_pois data is not covered by given interval ")
    return(NULL)
  }

  hist.pre <- hist(data,breaks = seq(interval[1],interval[2],length.out = n.split + 1 ),plot=FALSE)
  
  hist.res <- hist(hist.pre$counts,breaks=seq(-0.5,max(hist.pre$counts)+0.5,by=1),plot=FALSE)

  lambda = length(data)/n.split
  ch.res <- chisqr_genral(hist.res, function(x) ppois(x,lambda=lambda),doplot )

  return( list( chisqr= ch.res$chisqr , dof = ch.res$n.intervals ) )
}

chisqr_chisqr <- function(data,dof,doplot=FALSE){
  hist.res <- hist(data,plot=FALSE)
  ch.res <- chisqr_genral(hist.res, function(x) pchisq(x,df=dof),doplot )

  return( list( chisqr= ch.res$chisqr , dof = ch.res$n.intervals ) )
}


chisqr_check <- function(data,f,...,conf=0.95,doplot=FALSE){
  res <- f(data,...,doplot=doplot)
  return( pchisq(res$chisqr,df=res$dof) < conf )
}
