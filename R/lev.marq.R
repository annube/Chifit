


lev.marq <- function( x,y,dy,fn,dfn,par.0,par.priors=par.0,priors.sigmas=rep(Inf,length(par.0)),upd.range,maxit=100,verbose=F ){


  beta = par.0
  beta.new = beta
  lambda = 1.
  tol = 1.e-10
  nu = 2

  ## weight matrix
  W <- diag(1/dy^2)
  
  res = y - fn(x,beta)
  chisqr = sum( ( res/dy )^2 ) + sum((beta - par.priors)^2/priors.sigmas^2)


  chisqr.fn <- function(par) sum( ((y - fn(x,par)) /dy )^2 ) + sum((par - par.priors)^2/priors.sigmas^2)
  
  rel.diff = 1


  num.it = 0
  while( rel.diff > tol && num.it<=maxit){
    num.it = num.it + 1
    
    X <- dfn(x,beta)[,upd.range]
    rhs = t(X) %*% W %*% res - (beta-par.priors)/priors.sigmas^2
    Cw <- t(X) %*% W %*% X + diag( 1/priors.sigmas^2)
   
    if(verbose)    print( sprintf(" lambda = %e " ,lambda) )
    dbeta =  ( as.vector( lev.marq.solution(rhs,Cw,lambda) ) )

    ## primitive line search
##     obj.fn <- function(alpha) chisqr.fn(beta+alpha*dbeta)
##     or <- optim(c(0),obj.fn)

    
    
    beta.new[upd.range] = beta[upd.range] + dbeta
##    print(beta.new)
    res.new = try( y - fn(x,beta.new) )
    if( inherits( res.new , "try-error" ) || any( is.nan( res.new ) ) ) {
      print("increasing lambda because of bad evaluation of new residue " )
      lambda = lambda * nu
      next
    }
##    print(res.new)
    chisqr.new = sum( ( res.new/dy )^2 ) + sum((beta.new - par.priors)^2/priors.sigmas^2)

    if( chisqr.new < chisqr  ) {
      rel.diff = abs( 1-chisqr.new/chisqr )
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
        print(paste("increasing lambda because of increase of new chisqr (new chisqr = " , chisqr.new ,")" ))
        lambda = lambda * nu
    }
  }


  return( list( beta = beta , Chisqr = chisqr ) )
  
  
}


lev.marq.simple <- function( x,y,dy,fn,dfn,par.0,upd.range,maxit=100,verbose=F,lambda ){


  beta = par.0
  beta.new = beta
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
   
    dbeta =  ( as.vector( lev.marq.solution(rhs,Cw,lambda) ) )
    beta.new[upd.range] = beta[upd.range] + dbeta

    res.new = try( y - fn(x,beta.new) )
    
    chisqr.new = sum( ( res.new/dy )^2 )

    rel.diff = abs( 1-chisqr.new/chisqr )
    if(verbose){
      print( rel.diff )
      print( chisqr.new )
      print(beta.new)
    }
    beta = beta.new
    res = res.new
    chisqr = chisqr.new

  }


  return( list( beta = beta , Chisqr = chisqr ) )
  
  
}



lev.marq.solution <- function( rhs , Cw , lambda ) {

  Scale <- diag( 1/sqrt(diag(Cw) ))
  Cw.s <- Scale %*% Cw %*% Scale
  rhs.s <- as.vector( Scale %*% rhs )

  res <- try(
             solve( Cw.s + diag ( rep(lambda,length(rhs) ) ) ,rhs.s)
             )

  if( inherits( res , "try-error") )
    return( lev.marq.solution.svd( rhs.s,Cw.s,lambda ) )

  return( as.vector( Scale %*% res) ) 
}

lev.marq.solution.svd <- function(rhs,Cw,lambda) {
  svd.res <-  svd( Cw + diag ( lambda * diag(Cw ) ) )
  di = svd.res$d
  non.zero.range = ! di == 0 
  di[non.zero.range] = 1/di[non.zero.range]

  svd.res$v %*% ( diag(di) %*% ( t(svd.res$u) %*% rhs ) )
  
}

  
