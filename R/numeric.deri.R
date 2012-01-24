
## in case we dont want to calculate the derivative
numeric.deri <- function(f,x,tol=1.e-10){

  n <- length(x)
  df <- numeric(n)
  for( i in 1:n){

    eps = 1.
    delta=numeric(n)
    delta[i]=0.1

  
    df[i] <- (f(x+delta)-f(x))/delta[i]
    while(df[i] == 0 ){
      delta[i] = delta[i] * 10
      df[i] <- (f(x+delta)-f(x))/delta[i]
    }
    
    while( eps > (tol) ) {
      delta <- delta*0.5
      df_new <- (f(x+delta)-f(x))/delta[i]
      eps <- abs((df_new-df[i]))
      df[i] = df_new
      print(df)
    }
    print("=============")
  }
  return(df)
}


