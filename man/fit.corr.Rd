

\name{mcapply}
\alias{mcapply}
\docType{package}
\title{
multicore version of apply
}
\description{
  This function ...
}
\usage{
mcapply( X , fn , ...)
}

\arguments{
  \item{X}{A matrix with rows of data. Each row is applied to the function.}
  \item{fn}{The function to which the data rows should be applied to.}
}
\seealso{\code{\link{multicore::mclapply}}}
\examples{
  X <- array( rnorm(100) , dim=c(10,10) )
  fn <- function( y ) { x = 1:10 ; lm(y~x)$coefficients }
  X.lm.fit <- mcapply( X , fn )
}
