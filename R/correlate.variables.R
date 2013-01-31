

correlate.variables <- function( C , C.org=diag(rep(1,dim(C)[1])) ) {


  C.org.eig <- eigen(C.org)
  C.org.sqrt.i <- C.org.eig$vectors %*% diag( 1/sqrt(C.org.eig$values) )

  C.eig <- eigen(C)
  C.eig.sqrt <- diag(sqrt(C.eig$values) ) %*% t( C.eig$vectors )

  return( C.org.sqrt.i %*% C.eig.sqrt )
  
}
