bisect <- function(f,x1,x2,tol=1.e-7){

  y1 <- f(x1)
  y2 <- f(x2)

  m = (y2-y1)/(x2-x1)
  n = y2-m*x2
  
  
  mid <- - n/m

  fm <- f(mid)

  if(fm == 0 ) return(mid)
  
  if( abs(fm) / min(c(abs(f(x1)),abs(f(x2))) ) > 0.99 )
    return(mid)
  
  if( fm > 0 )
    return (bisect(f,x1,mid))
  else
    return (bisect(f,mid,x2))
  
}


ridder <- function(f,x1,x2,tol=1.e-7,maxit=100){
  f1 <- f(x1)
  f2 <- f(x2)
  

  if( ( f1 < 0 && f2 > 0 ) ){
    res <- Inf
    for( i in 1:maxit ){
#      print(paste(" i " , i))
        mid <- 0.5*(x1+x2)
        fmid <- f(mid)

        s <- sqrt(fmid^2-f1*f2)
        if(s == 0 ) return(res)
        xnew <- mid+(mid-x1)*sign(f1-f2)*fmid/s
        
        if( abs(xnew-res) <= tol ) return(res)
        res <- xnew
        fnew=f(res)
        if(fnew == 0) return(res)

        if(fnew > 0 ){
          f2 <- fnew
          x2 <- res
        } else {
          f1 <- fnew
          x1 <- res
        }
        if( abs(x1-x2) <= tol ) return (res)
    }
    print("maximum number of iterations reached")
      
    
  } else {
    if(f1==0) return (x1)
    if(f2==0) return (x2)
    print("Error root must be bracketed")
  }
}
