


newton <- function(f, df, x0, ..., maxit=100,tol=1.e-10) {

  x=x0
##  y = f(x,...)

  dx = 0
  first = TRUE
  it = 0
  while( ( abs(dx/x) > tol && it < maxit ) || first ) {
    first = FALSE
    y = f(x,...)
    dx=-y/df(x,...)
    x=x+dx
##    print( dx/x)
    if( is.nan( dx/x ) ){
      x = rnorm(1)
      dx = 0
      first = TRUE
    }
    it = it + 1
  }
##  if( it == maxit ) print("warning maxit reached ")
  return(x) 
}
