



block.gen.test.data <- function(n,delta=0.1){
  x0 = -2

  data = numeric(n)

  x = x0
  
  for( i in 1:n ){
    newx = x + delta * ( 2. * runif(1) -1. )
    dS = (newx^2)/2 - (x^2)/2
    if( exp(-dS) > runif(1) ) {
      x = newx
    }
    data[i] = x
  }

  return(data)
  
}


block.to.half <- function(x) {
  half = floor( length(x)/2 )
  if( half > 0 )
    return( 0.5 * ( x [ 2 * (1:half) - 1 ] + x [ 2 * (1:half) ] ) )
  else
    return( numeric(0) )
}


log2.trunc <- function( x ) {
  e = 0
  while( x > 1 ){
    x = floor( x / 2 )
    e = e + 1
  }
  return ( e )
}



blocking.analysis <- function(x) {

  sigma.x = numeric( log2.trunc(length(x)) - 1  )
  dsigma.x = numeric( log2.trunc(length(x)) - 1 )

  for( i in 1:( log2.trunc( length(x) ) -1 ) ) {
    sigma.x[i] = sqrt( var( x )/length(x) )
    dsigma.x[i] = sigma.x[i] / sqrt(2 * ( length(x) - 1 ) )
    x <- block.to.half( x )
  }
  
  return( cbind( sigma.x , dsigma.x ) )
}



calc.gamma.t <- function(x,t) {
  n = length(x)
  x.mean=mean(x)
  return( mean( ( x[1:(n - t )] - x.mean  ) * ( x[(t+1):n] - x.mean ) ) )
}

correl.gamma.t <- function(x,T = 10 ) {
  return( sapply(0:T , function(t) calc.gamma.t(x,t) ) )
}
  

error.gamma.t <- function(data,T) {
  gamma.t <- correl.gamma.t ( data,T)
  ( gamma.t[1] + 2 * sum(  (1 - (1:T)/length(data)) * gamma.t[2:(T+1)] ) ) /
    ( length(data) - 2 * T - 1 + T*( T + 1) / length(data) )
}
  
