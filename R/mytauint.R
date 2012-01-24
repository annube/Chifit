

mytauint <- function(x,gammalen=10) {
  
  len <- length(x)
  numwindows=10
  wl <- floor( (len-gammalen)/numwindows )

  

  Gamma <- array(0,dim=c(gammalen,numwindows))

  for ( window in 0:(numwindows-1) )
    for( i in 0:(gammalen-1) ) {
      for( j in 0:(wl-1) ){
        xindex = 1+wl*window+j
        Gamma[i+1,window+1]=Gamma[i+1,window+1]+x[xindex]*x[xindex+i]
      }
  }

  

  return(Gamma)
  
}
