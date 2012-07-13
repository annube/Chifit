


make_hr_quantity <- function(q,dq,min_err=5,num_digits) {


  if( is.na(q) || is.na(dq) )
    return(sprintf("NA(NA)"))
  
  if( missing(num_digits) && ! is.na(dq) && ! is.null(dq) ) {
    if( abs(dq) < 1.e-6) {
      num_digits = 2
    } else {
      num_digits = ceiling( log(min_err/dq)/log(10) )
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
