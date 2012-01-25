



# load hadron library

library(hadron)


## partitions <- read.table("partitions.csv",header=TRUE)


data(partitions)


g.tilde<- function(lambda,partitions,maxN=20)  {
  
  result=0
  for( i in 1:maxN){
    x = sqrt(i)*lambda
    index=which(partitions$n == i)
    m = partitions$m[index]
    result=result+4*m/x*besselK(x,1)
  }
  return(result)
}



g.tilde.deri<- function(lambda,partitions,maxN=20)  {
  
  result=0
  for( i in 1:maxN){
    x = sqrt(i)*lambda
    index=which(partitions$n == i)
    m = partitions$m[index]
    result = result +
     4 * m * (  -1. / ( x * lambda ) * besselK( x, 1 ) - 0.5 * sqrt( i ) / x * ( besselK( x, 0 ) + besselK( x, 2 ) ) )
  }
  return(result)
}


g.tilde.lookup.len = 2000
g.tilde.lookup = array(dim=c(g.tilde.lookup.len,3))

g.tilde.lookup.max=10
g.tilde.lookup.min=1.e-10

g.tilde.lookup[,1] <- seq( g.tilde.lookup.min , g.tilde.lookup.max, length.out=g.tilde.lookup.len)
g.tilde.lookup[,2] <- g.tilde(g.tilde.lookup[,1],partitions)
g.tilde.lookup[,3] <- g.tilde.deri(g.tilde.lookup[,1],partitions)



g.tilde_<- function(lambda,partitions,maxN=20)  {

  result <- numeric(length(lambda))
  
  for( i in 1:length(lambda) ){
    if( lambda[i] > g.tilde.lookup[1,1] && lambda[i] < g.tilde.lookup[g.tilde.lookup.len,1] ){
##      print(paste("looking up .."))#
      absdiff.lu <- abs(g.tilde.lookup[,1]-lambda[i])
      indices <- order(absdiff.lu)[1:3]

      if( lambda[i] == g.tilde.lookup[indices[1],1])
        result[i] = g.tilde.lookup[indices[1],2]
      else  if( lambda[i] == g.tilde.lookup[indices[2],1])
        result[i] = g.tilde.lookup[indices[2],2]
      else {
##         alpha <- (g.tilde.lookup[indices[1],1] - lambda[i]) / (g.tilde.lookup[indices[1],1] - g.tilde.lookup[indices[2],1])
##         result[i] <- (1-alpha) * g.tilde.lookup[indices[1],2] + alpha * g.tilde.lookup[indices[2],2]

        M <- cbind( g.tilde.lookup[indices[1:3],1]^2 , g.tilde.lookup[indices[1:3],1] , rep(1,3) )
        pl.co <- solve(M,g.tilde.lookup[indices[1:3],2])
        result[i] <- sum( pl.co * c(lambda[i]^2,lambda[i],1) )
      }

      
    } else {
      print("not looking up in g.tilde")
      result[i] <- g.tilde(lambda[i],partitions,maxN)
    }
  }
  return( result )
}

g.tilde.deri_<- function(lambda,partitions,maxN=20)  {

  result <- numeric(length(lambda))
  
  for( i in 1:length(lambda) ){
    if( lambda[i] > g.tilde.lookup[1,1] && lambda[i] < g.tilde.lookup[g.tilde.lookup.len,1] ){
##      print(paste("looking up .."))#
      absdiff.lu <- abs(g.tilde.lookup[,1]-lambda[i])
      indices <- order(absdiff.lu)[1:3]

      if( lambda[i] == g.tilde.lookup[indices[1],1])
        result[i] = g.tilde.lookup[indices[1],3]
      else  if( lambda[i] == g.tilde.lookup[indices[2],1])
        result[i] = g.tilde.lookup[indices[2],3]
      else {
##         alpha <- (g.tilde.lookup[indices[1],1] - lambda[i]) / (g.tilde.lookup[indices[1],1] - g.tilde.lookup[indices[2],1])
##         result[i] <- (1-alpha) * g.tilde.lookup[indices[1],3] + alpha * g.tilde.lookup[indices[2],3]

        M <- cbind( g.tilde.lookup[indices[1:3],1]^2 , g.tilde.lookup[indices[1:3],1] , rep(1,3) )
        pl.co <- solve(M,g.tilde.lookup[indices[1:3],3])
        result[i] <- sum( pl.co * c(lambda[i]^2,lambda[i],1) )
      }

      
    } else {
      print("not looking up in g.tilde")
      result[i] <- g.tilde.deri(lambda[i],partitions,maxN)
    }
  }
  return( result )
}

