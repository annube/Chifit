


make_hr_quantity <- function(q,dq,min_err=5,num_digits) {


  if( is.na(q) || is.na(dq) )
    return(sprintf("NA(NA)"))
  
  if( missing(num_digits) && ! is.na(dq) && ! is.null(dq) ) {
    if( abs(dq) < 1.e-6) {
      num_digits = 2
    } else {
      num_digits = max(0,ceiling( log(min_err/dq)/log(10) ))
    }
  } else {
    if( missing( num_digits ) ) {
      num_digits = 1

      q = 0
      dq = 0
    }
  }
  

  format=sprintf("%%%d.%df(%%d)",(num_digits+3),num_digits)

  return( sprintf(format,q,round(dq*10^num_digits)) )
}


## the inverse function

convert.from.hr.quantity <- function(string,expected.len=1) {
  ## remove whitespace
  string <-  gsub(" " , "" , string)

  ## split away different parts
  parts <- unlist( strsplit(string,"[()]") )

  non.zero.len.range <- which( sapply(parts,nchar) != 0 )

  parts = parts[non.zero.len.range]

  
  
  print( parts )
  ## determine number of digits after the comma

  after.comma.string <- unlist ( strsplit( parts[1] , split = "[.]" ) ) [2]
  digits <- length( unlist( strsplit(after.comma.string,split="") ) )

  value <- try( as.numeric( parts[1] ) ,TRUE)

  res <- numeric(length(parts))
  res[1] = value

  if( length(parts) > 1 )
    for( i in 2:length(parts) ){
      digits.err <- length( unlist( strsplit(parts[i],split="") ) )
      error.str <- paste0( c( "0." , rep("0" , digits-digits.err ) , parts[i]) , collapse="" )
      print(error.str)
      res[i] = as.numeric(error.str)
    }

  if( is.na( res[1] ) )
    return( rep( NA ,expected.len ) )
  else 
    return( res )
  
}                             
