library(chifit)
expfn<-function(x,par,...) par[1] * exp(-par[2]*x)
dexpfn<-function(x,par,...) cbind(  exp(-par[2]*x) , -par[1]*x*exp(-par[2]*x) )


coshfn<-function(x,par,aargs) par[1] * cosh( par[2] * (x - aargs$Thalf ) )
dcoshfn<-function(x,par,aargs) cbind(
                                     cosh( par[2] * (x - aargs$Thalf ) ) ,
                                     par[1] * ( x - aargs$Thalf ) * sinh( par[2] * (x - aargs$Thalf ) )
                                     )



 tab.res <- read.table("obs_s2_N200_a0.050000_M00.500000_musq0.000000_l1.000000_J0.000000.csv",header=TRUE,colClasses = "numeric",nrows = 100000 )


Tmax = length(tab.res) - 6

Corr <- as.matrix( tab.res[,0:Tmax + 6] )

avg.corr.ts <- function(Corr) apply(Corr,2,mean)


boot.R = 1000

boot.res <- tsboot( Corr , avg.corr.ts, boot.R, l = 10 , sim="fixed" ,parallel = "multicore",ncpus = 16)

C <- cov( boot.res$t)


range = 15:41


make.fit <- function( Corr ) {
  Corr.avg <- apply( Corr, 2, mean )
  return(
   wnlls(x = 0:Tmax ,
         y = Corr.avg,
         ##dy = sqrt( diag( C ) ),
          C = C,
         f= coshfn,
         df = dcoshfn,
         par = c(0.1,0.1),
         range = range,
         aargs = list( Thalf = 100 )
         )
         )
}

wnlls.res <- make.fit( Corr )
plot.wnlls(wnlls.res)


make.fit.tsboot <- function( Corr ) {
  wnlls.res <- make.fit( Corr )
  return( c( wnlls.res$beta, wnlls.res$Chisqr ) )
}

make.fit.boot <- function( Corr,i ) {
  wnlls.res <- make.fit( Corr[i,] )
  return( c( wnlls.res$beta, wnlls.res$Chisqr ) )
}

boot.R.fit <- 1000
##boot.res.fit <- tsboot( Corr , make.fit.boot, boot.R.fit, l = 10 , sim="fixed" ,parallel = "multicore",ncpus = 16)
boot.res.fit <- boot( Corr , make.fit.boot, boot.R.fit ,stype="i", parallel = "multicore",ncpus = 14)
