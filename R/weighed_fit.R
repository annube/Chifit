

# load hadron library

library(hadron)


truncate.mat <- function(C,numev =dim(C)[1]){
  C.eig <- eigen(C)
  dD <- rep(0,dim(C)[1])
  dD[1:numev] =  C.eig$values[1:numev]
  return( (C.eig$vectors) %*% diag(dD) %*% t(C.eig$vectors) )
}


## weighed least squares fit with gaussian error calculation (propagation)
wlm <- function(x,y,dy=NULL,C=NULL){

  ## error_propagation weighed

  

  if( ! is.null(dy) ){  
    w <- 1/dy^2
    W <- diag(w)
  }

  if ( ! is.null(C) ) {
    W = solve(C)
    if( is.null(dy) ) dy <- sqrt(diag(C))
  }
  if( is.null(dy) && is.null(C) ) {
    print("Error in wlm if dy is not given, give at least the correlation matrix of the y's")
    return(NULL)
  }
  
  if( ! is.null(x) ) {
    
    if( ! is.null( dim (x) ) )
      len <- length(x[,1])
    else
      len <- length(x)
    
    X <- cbind(  x )
#    X <- cbind( rep(1,len) , x )
    
  }
  else {
    X <- cbind( rep(1,length(y)) )
  }

#  print(dim(X))

  dof=len-dim(X)[2]

#  print(X)
#  print(W)
  
  Cw <- t(X) %*% W %*% X

  
  Cwi <- solve(Cw)

  
  Beta <-  (  Cwi %*%  (  t(X) %*% ( W  %*% y ) ) )
  
## this formula comes from the error propagation
  dBetaw <- sqrt(  ( Cwi %*%  t(X) %*% W )^2  %*% (dy)^2  )



  yh <-  X %*% Beta
  delta.y=as.matrix(y)-yh
  chisqr <- sum ( ( W %*% delta.y ) * delta.y )


  return(list(dbeta=dBetaw,beta=Beta,
              x=as.matrix(x),y=y,dy=dy,
              predict = yh,
              dof=dof,
              dpredict = sqrt( (X^2) %*% (dBetaw^2) ),
              Chisqr=chisqr
              )
         )
}


expfn <- function(x,par){
  par[1]*exp(par[2]*x)
}

exp2fn <- function(x,par){
  par[1]*exp(par[2]*x) + par[3]*exp(par[4]*x)
}

dexpfn <- function(x,par){
  cbind(exp(par[2]*x) , x*par[1]*exp(par[2]*x) )
}

dexp2fn <- function(x,par){
  cbind(exp(par[2]*x) , x*par[1]*exp(par[2]*x) , exp(par[4]*x), x*par[3]*exp(par[4]*x) )
}


chisqrfn <- function(par,x,y,fnctn,dfnctn,W,...){
##  print(par)
  r <- as.vector( y-fnctn(x,par,...) )
  res <- sum( r * ( W %*% r ) )
##  print(res)
  return(res)
}

dchisqrfn <- function(par,x,y,fnctn,dfnctn,W,...){
  r <- as.vector( y-fnctn(x,par,...) )
  df <-  dfnctn(x,par,...)
  apply( df , 2 , function(x) -2*sum( x * ( W %*% r ) ) )
}


## weighed least squares fit with gaussian error calculation (propagation)
wnlls <- function(x.in,y,dy=NULL,C=NULL,f,df=NULL,par,range=NULL,method="optim",optim.method="BFGS",aargs=NULL,...){

  ## error_propagation weighed


  x=as.matrix(x.in)

  if(is.null(range) )
    range=1:dim(x)[1]

  if( ! is.null(dy) ){  
    w <- 1/dy[range]^2
    W <- diag(w)
  }

  if ( ! is.null(C) ) {
##    print("trying to solve C")
    print(C[range,range])
    W = solve(C[range,range])
    if( is.null(dy) )
      dy = sqrt(diag(C))
  }
  if( is.null(C) && is.null(dy) ) {
    print("Error in wlm if dy is not given, give at least the correlation matrix of the y's")
    return(NULL)
  }




  if(is.null(df)) {
    method="optim"
##    library(maxLik)
##    df=function(x,par,aargs) numericGradient(function(x_)  f(x,x_,aargs) ,t0=par,eps=1.e-10  )
  }

## do solution
  if( method == "mycg" && !is.null(df) ){
    print("starting mycg optimisation")
    optim.res <- mycg(par,x[range,],y[range],f,df,W,tol=1.e-10,maxsteps=100,...)
    print("finished: mycg optimisation")
    beta=optim.res$beta
    num.it=optim.res$num.it
    betahist=optim.res$betahist
  } else {
    
    if(is.null(df))
      grad=NULL
    else
      grad=dchisqrfn
    
    optim.res <- optim(par,chisqrfn,
                       gr=grad,
                       x=x[range,],y=y[range],W=W,fnctn=f,dfnctn=df,aargs=aargs,...,method=optim.method
                       )



    beta = optim.res$par
##    print(beta)
    num.it = 1
    betahist=beta
  }

  
  delta.y=as.vector( y[range] - f(x[range,],beta,aargs=aargs) )
  chisqr <- sum ( ( W %*% delta.y ) * delta.y )

  ## calculate prediction of the model and error of it
  predict=as.vector( f(x,beta,aargs=aargs) )

  dof <- length(range) - length(beta)
  fit.consistend.chisq <- ( chisqr < qchisq(0.95,df=dof) )

  if( is.null(df) ){
    res <- list(
                ## calling parameters
                x=x,y=y,dy=dy,C=C,f=f,par=par,range=range,
                ##
                beta=beta,
                predict=predict,
                dpredict=rep(0,length(predict)),
                ##
                Chisqr = chisqr,
                dof = dof,
                fit.consistend.chisq = fit.consistend.chisq
                )
    
    
    res$optim.res = optim.res
    
    return(res)
  }
  
  J=df(x[range,],beta,aargs=aargs)


  Cw <- t(J) %*% W %*% J
  
  Cwi <-  try( qr.solve(Cw,tol=1.e-20) )
  while(inherits(Cwi,"try-error") ) {
    beta = beta+0.1*rnorm(length(beta))
    J=df(x[range,],beta,aargs=aargs)
    Cw <- t(J) %*% W %*% J
    Cwi <-  try( qr.solve(Cw,tol=1.e-20) )
  }
  
   M <- Cwi %*% t(J) %*% W

  dbeta <- as.vector( sqrt( M^2 %*% dy[range]^2 ) )
  

  Jf=df(x,beta,aargs)


  ## calculate prediction of the model and error of it
  dpredict=as.vector( sqrt( (Jf^2) %*% (dbeta^2) ) )
  
  ## check consistency of predicted and observed data within error bounds
  
  norm.diff <- (predict-y)/dpredict
  crit.fullfilled <- sum( abs(norm.diff[range])<1.96 )
  fit.consistend.pred = ( crit.fullfilled == length( norm.diff[range] ) )

  ## investigate linearity of the function in the fitted Region

  linearity.beta=logical(length(beta))
  for( betai in 1:length(beta)){

    incbeta <- numeric(length(beta))
    incbeta[betai] <- 2*dbeta[betai]
    
    rel.lin.approx.error.pos <- as.vector(
                                          f(x[range,],beta+incbeta,aargs=aargs)
                                          - ( f(x[range,],beta,aargs=aargs) + ( J %*% incbeta ) )
                                          ) /y[range]
    rel.lin.approx.error.neg <- as.vector(
                                          f(x[range,],beta-incbeta,aargs=aargs)
                                          - ( f(x[range,],beta,aargs=aargs) - ( J %*% incbeta ) )
                                          ) /y[range]

    ## we allow violations of the linearity up to 1%
    linearity.beta[betai] <- ( max ( abs(rel.lin.approx.error.pos) ) < 0.01  &&
                       max ( abs(rel.lin.approx.error.neg) ) < 0.01 )
    
  }
  
  
  
  res <- list(
               ## calling parameters
               x=x,y=y,dy=dy,C=C,f=f,df=df,par=par,range=range,
               ##
               beta=beta,dbeta=dbeta,betahist = betahist,
               predict=f(x,beta,aargs=aargs),
               dpredict=sqrt( (Jf^2) %*% (dbeta^2) ),
               ##
               num.it = num.it,
               Chisqr = chisqr,
               dof = dof,
               fit.consistend.chisq = fit.consistend.chisq,
               fit.consistend.pred = fit.consistend.pred,
              linearity.beta=linearity.beta

               )

  
  res$optim.res = optim.res

  return(res)


}


plot.wnlls <- function(res,use.col=1,x.data = NULL,plot.range,...){


##  attach(res)

  if( is.null(res$dy) ) res$dy = sqrt(diag(res$C) )

  if( is.matrix(res$x) )
    x= res$x[,use.col]
  else
    x = res$x

  if( !missing(x.data) ){
    x = x.data
  }

  if(missing(plot.range))
    plot.range = res$range
  
  plotwitherror(x[plot.range],res$y[plot.range],res$dy[plot.range],...)

  
  plotwitherror(x[res$range],res$y[res$range],res$dy[res$range],rep=TRUE,col="green")

  
  
  title(main=bquote(paste("Fit of " , y ) ))

  x.min <- min(x[plot.range])
  x.max <- max(x[plot.range])

  lines( x[plot.range],res$predict[plot.range] )
  lines( x[plot.range],(res$predict+res$dpredict)[plot.range],lty="dashed" )
  lines( x[plot.range],(res$predict-res$dpredict)[plot.range],lty="dashed" )
##  detach(res)



  text(  0.5*(x.min+x.max), min( ( res$y - res$dy ) [plot.range] ),  bquote(  paste ( chi^2 == .( res$Chisqr ) , "; " , dof == .(res$dof) ) ) ) 
  text(  0.5*(x.min+x.max), max( ( res$y + res$dy ) [plot.range] ),  bquote(  paste ( beta[1] == .( res$beta[1] ) +- .( res$dbeta[1] ),
                                                                "; " ,
                                                                beta[2] == .(res$beta[2]) +- .( res$dbeta[2] )) ) )

  
}

wnlls.test1 <- function(par=c(1,0),...){

  corr <- read.table("R/corr_ana_N80_a0.250000_M01.000000_cos_s1.csv",header=TRUE)
  corr.cov <- as.matrix(
                        read.table("R/corr_cov_N80_a0.250000_M01.000000_cos_s1.csv",header=TRUE)
                        )

  skip.1st <- 4

  range=(skip.1st+1):length(corr$T)
  
##  res <- wnlls(corr$T,corr$Corr,C=corr.cov,f=expfn,df=dexpfn,par=c(1.0,0.0),range=range)

  res <- wnlls(corr$T,corr$Corr,dy=corr$DCorr,f=expfn,df=dexpfn,par=par,range=range,...)


  
  plot.wnlls(res)
  
  return(res)
}


wnlls.test2 <- function(...){

  corr <- read.table("R/corr_ana_N80_a0.250000_M01.000000_cos_s1.csv",header=TRUE)
  corr.cov <- as.matrix(
                        read.table("R/corr_cov_N80_a0.250000_M01.000000_cos_s1.csv",header=TRUE)
                        )

  skip.1st <- 7

  range=(skip.1st+1):length(corr$T)
  

  xy <- make_grid(x=seq(0,5,length.out=10),y=seq(-1.,0,length.out=10) )

  xy <- cbind(xy$x,xy$y)
  
  
  applyfn <- function(x) {
    print(x)
    wnlls(corr$T,corr$Corr,dy=corr$DCorr,f=expfn,df=dexpfn,par=x,range=range,...)$beta
  }

  betas <- apply(xy,1,applyfn)

  res1 <- wnlls.test1(par=c(1,-1),method="mycg")
  

  
  return(list(xy=xy,betas=betas))
}


wlm.test <- function(){
  m = 2
  n = 1
  N = 10
  x <-  sort(rnorm(N))
  dy = 0.5 * abs(rnorm(N))
  y <-  m * x + n + dy*rnorm(N)

  wlm.res <- wlm(cbind(  x, rep(1,N)) ,y,dy)
  plot.wnlls( wlm.res )
}


df.regression.fit <- function(xs,ys,dy,f,df,par,n,...){
  for( k in 1:n) {
    res = ys - f(xs,par,...)
    dfv <- df(xs,par,...)
    print(par)
    lm.res = lm(res~dfv+0,weights = 1/dy^2)
    par = par + lm.res$coefficients
  }
  return(par)
}
