
\name{calc.deri.par.fn}
\alias{calc.deri.par.fn}
\docType{package}
\title{
Title
}
\description{
  This function calculates the derivative of a expression containing par
  and x as variables with respect to the members of the array par.
}
\usage{
calc.deri.par.fn( expr, parlen )
}
\examples{
f <- function(par,x) par[1]*exp(-par[2] * x )
df.expr <- calc.deri.par.fn( body(f)[[1]] , 2)
df <- ....
}
