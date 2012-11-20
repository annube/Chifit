


lev.marq <- function( x,y,dy,fn,dfn,par.0,upd.range,maxit=100,verbose=F ){


  beta = par.0
  beta.new = beta
  lambda = 1.
  tol = 1.e-10
  nu = 2

  ## weight matrix
  W <- diag(1/dy^2)
  
  res = y - fn(x,beta)
  chisqr = sum( ( res/dy )^2 )
  rel.diff = 1


  num.it = 0
  while( rel.diff > tol && num.it<=maxit){
    num.it = num.it + 1
    
    X <- dfn(x,beta)[,upd.range]
    rhs = t(X) %*% W %*% res
    Cw <- t(X) %*% W %*% X
    if(verbose)    print( sprintf(" lambda = %e " ,lambda) )
    dbeta =  ( as.vector( lev.marq.solution(rhs,Cw,lambda) ) )
    beta.new[upd.range] = beta[upd.range] + dbeta
##    print(beta.new)
    res.new = try( y - fn(x,beta.new) )
    if( inherits( res.new , "try-error" ) || any( is.nan( res.new ) ) ) {
      lambda = lambda * nu
      next
    }
##    print(res.new)
    chisqr.new = sum( ( res.new/dy )^2 )

    if( chisqr.new < chisqr ) {
      rel.diff = 1-chisqr.new/chisqr
      if(verbose){
        print( rel.diff )
        print( chisqr.new )
        print(beta.new)
      }
      beta = beta.new
      res = res.new
      chisqr = chisqr.new
      lambda = lambda / nu
    } else {
      lambda = lambda * nu
    }
  }


  return( list( beta = beta , Chisqr = chisqr ) )
  
  
}


lev.marq.solution <- function( rhs , Cw , lambda ) {
  res <- try(
             solve( Cw + diag ( lambda * diag(Cw ) ) ,rhs)
             )

  if( inherits( res , "try-error") )
    return( lev.marq.solution.svd( rhs,Cw,lambda ) )

  return( res) 
}

lev.marq.solution.svd <- function(rhs,Cw,lambda) {
  svd.res <-  svd( Cw + diag ( lambda * diag(Cw ) ) )
  di = svd.res$d
  non.zero.range = ! di == 0 
  di[non.zero.range] = 1/di[non.zero.range]

  svd.res$v %*% ( diag(di) %*% ( t(svd.res$u) %*% rhs ) )
  
}

  
