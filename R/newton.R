


newton <- function(f, df, x0, ..., maxit=100,tol=1.e-10) {

  x=x0
  y = f(x,...)

  dx = 0
  first = TRUE
  it = 0
  while( ( abs(dx/x) > tol && it < maxit ) || first ) {
    first = FALSE
    dx=-y/df(x,...)
    x=x+dx
    y = f(x,...)
##    print( dx/x)
    it = it + 1
  }
  return(x) 
}
