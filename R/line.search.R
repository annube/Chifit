

x0.quad.regr <- function( x, fx ){
  (
   ( x[2]^2 - x[3] ^ 2 ) * fx[1] +
   ( x[3]^2 - x[1] ^ 2 ) * fx[2] +
   ( x[1]^2 - x[2] ^ 2 ) * fx[3]
   ) / (
        2 * (
             ( x[2] - x[3] ) * fx[1] +
             ( x[3] - x[1] ) * fx[2] +
             ( x[1] - x[2] ) * fx[3]
             )
        )
}


line.search <- function( f , x0 ,delta0,tol=1.e-10,maxit=10){


  x <- c(0,0,0)
  fx <- c(0,0,0)
  
  delta <- delta0
  x[2] <- x0
  fx[2] <- f(x[2])
  
  x[3] <- x0+delta
  fx[3] <- f(x[3])

  n <- 1
  while( fx[3] < fx[2] ){
    x[3]= x0 + delta * 2^n
    fx[3] <- f(x[3])
    n=n+1
  }


  x[1] <- x0-delta
  fx[1] <- f(x[1])
  n <- 1
  while( fx[1] < fx[2] ){
    x[1]= x0 - delta * 2^n
    fx[1] <- f(x[1])
    n=n+1
  }

  

##  print( paste( x) )
##  print( paste( fx ) )

  it <- 1
  while( x[3]-x[1]>tol && it < maxit ){
    it = it + 1
    new.x2 = x0.quad.regr(x,fx)


    nfx2 <- f(new.x2)
    
    if(new.x2 < x[2] ) {
      if( nfx2 < fx[2] ){
         x[3] = x[2]
         fx[3] = fx[2]
         x[2] = new.x2
         fx[2] = nfx2
       } else {
         x[1] = new.x2
         fx[1] = nfx2
       }
    } else {
      if( nfx2 < fx[2] ){
        x[1] = x[2]
        fx[1] = fx[2]
        x[2] = new.x2
        fx[2] = nfx2
      }else {
        x[3] = new.x2
        fx[3] = nfx2
      }
    }
##    print( paste( x[2],fx[2] ) )

##    print( fx[2] < fx[1] && fx[2] < fx[3] ) 
    
  }
  
  return(x[2])
}
