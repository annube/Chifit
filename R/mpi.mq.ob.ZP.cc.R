

mpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_mq_ob",x,par,aargs,F,PACKAGE="chifit")
}


dmpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_mq_ob_zp",x,par,aargs,T,PACKAGE="chifit")
}

## the "ZP"-versions of the following functions will be implemented in future
## mpi.0.mq.ob.cc <- function(x,par,aargs){
##   .Call("mpi_0_mq_ob",x,par,aargs,F,PACKAGE="chifit")
## }


## dmpi.0.mq.ob.cc <- function(x,par,x,aargs){
##   .Call("mpi_0_mq_ob",x,par,aargs,T,PACKAGE="chifit")
## }


## fpi.mq.ob.cc <- function(x,par,aargs){
##   .Call("fpi_mq_ob",x,par,aargs,F,PACKAGE="chifit")
## }


## dfpi.mq.ob.cc <- function(x,par,aargs){
##   .Call("fpi_mq_ob",x,par,aargs,T,PACKAGE="chifit")
## }


