# P_n(x) = params[1] x^n + ... + params[n-1] x + params [n]

## poly<-function(params,x_vals){
##   retval=rep(params[1],length(x_vals))
##   for(i in 2:length(params) ){
##    retval=retval*x_vals+params[i]
##   }
##   return(retval)
## }


## dpoly <- function(params,x_vals){

##   degree <- length(params) - 1

##   retval = numeric( degree + 1 )

##   retval = array( dim = c(degree + 1 , length( x_vals ) ) )

##   retval[(degree+1),] = 1

##   for( i  in (degree):1 )
##    retval[i,] = retval[(i+1),] * x_vals
##    return(retval)
## }

# y = sqrt( params[1]^2 - (x-parmas[2])^2 ) + params[3]

circle <- function(params,x_vals) {
  sqrt( params[1]^2 - (x_vals-params[2])^2 ) + params[3]
}


# derivative of circle

dcircle <- function(params,x_vals) {

  sqr = sqrt( params[1]^2 - (x_vals-params[2])^2 )

  retval=array(dim=c(3,length(x_vals)))
  retval[1,]=params[1]/sqr
  retval[2,] = (x_vals-params[2])/sqr
  retval[3,] = 1

  return( retval )
}

## derivative of circle

dcircledx <- function(params,x_vals) {

  sqr = sqrt( params[1]^2 - (x_vals-params[2])^2 )

  return( ( params[2] - x_vals ) / sqr )
}


# square root

square_root <- function(params,x_vals){

  params[1] * sqrt(x_vals-params[2]*0.01 ) + params[3]*0.01
}

# circle middlepoint determination
# need at least 3 points
calc_circle_params <- function(xs,ys){

 N=length(xs)

 mi=numeric( (N*(N-1))/2 )
 ni=numeric( (N*(N-1))/2 )

 for(i  in 1:(N-1) ){
  for(j  in (i+1):N ){

    index= j - i + N * ( N - 1 ) / 2  - ( N - i ) * ( N - i + 1 ) / 2

 #    print(paste(i,j,index))

    mi[index] = - 1 / ( ( ys[j] - ys[i] )  / ( xs[j]-xs[i] ) )
    ni[index] = 0.5 *  ( ys[j] + ys[i] ) - mi[index] * 0.5 * ( xs[j] + xs[i] )
  }
 }


 smi=sum(mi)
 smisq=sum(mi^2)
 sni=sum(ni)
 smini=sum(mi*ni)

 M= matrix(c(smisq , -smi, 
              -smi , N*(N-1)/2 ),
            nrow=2,ncol=2,byrow=TRUE)

 middle=solve(M)%*%c(-smini,sni)
 
# R =  max( sqrt( (middle[1]-xs)^2 + (middle[2]-ys)^2  ) ) 
 R =  max( sqrt( mean((middle[1]-xs)^2) + mean( (middle[2]-ys)^2 )  ) ) 

# test_res<-circle(c(R,middle),xs)
# indices_nan <- which( is.nan(test_res) )
# if( length( indices_nan ) > 0 )
#   R=sqrt(  

 return(c(R,middle))
 
}



chisqr_mk_function <- function ( param , input_vals) {

   sum( ( ( input_vals$fnctn ( param , input_vals$data_x ) - input_vals$data_y ) / input_vals $ data_dy ) ^ 2 )

}

dchisqr_mk_function <- function ( param , input_vals) {


   df = numeric(length(param))

  dims = length(dim(input_vals$data_y))

  dif = input_vals$dfnctn(param,input_vals$data_x)

  if(dims == 1) {
  for( i in 1:length(param) ){

    df[i] =  sum( 2 *  dif[i,] *
                 ( input_vals$fnctn ( param , input_vals$data_x ) - input_vals$data_y ) / input_vals $ data_dy^2   )
  }
  } else if(dims == 2) {
   for( i in 1:length(param) ){

    df[i] =  sum( 2 *  dif[i,,] *
                 ( input_vals$fnctn ( param , input_vals$data_x ) - input_vals$data_y ) / input_vals $ data_dy^2   )
  }
  } 


  return(df)

}


