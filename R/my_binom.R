
my_choose <- function(n,x) {
  e=exp(1)
  res = sqrt(n/(2*pi*x*(n-x) )) * (n/x)^x * (n/(n-x))^(n-x) *exp(1/12*(1/n-1/x-1/(n-x)) ) 
# 4.th order #*exp(-1/360*(1/n^3-1/x^3-1/(n-x)^3) )
  return(res) 
}

my_choose_ref <- function(n,x) {
  
  return( gamma(n+1)/gamma(x+1)/gamma(n-x+1) ) 
}


my_dbinom<-function(x,prob,size){
  my_choose_ref(size,x) * prob ^x * (1-prob)^(size-x)
}


