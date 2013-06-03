

jackknife <- function(x, blocklen = 10) {
  n = length(x)
  blocks = floor(n / blocklen)
  means = sapply( 0:(blocks-1) , function(i ) mean( x[1:blocklen + i * blocklen ] ) )
  M <- toeplitz( c(0,rep(1,blocks-1)))/(blocks-1)
  print(M)
  return( as.vector( M%*%   means ) )
  
}
