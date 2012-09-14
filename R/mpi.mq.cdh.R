

mpi.mq.cdh.cc <- function(x,par,aargs.par){
  .Call("mpi_mq_cdh",x,par,aargs.par,F,PACKAGE="chifit")
}

dmpi.mq.cdh.cc <- function(x,par,aargs.par){
  .Call("mpi_mq_cdh",x,par,aargs.par,T,PACKAGE="chifit")
}

fpi.mq.cdh.cc <- function(x,par,aargs.par){
  .Call("fpi_mq_cdh",x,par,aargs.par,F,PACKAGE="chifit")
}

dfpi.mq.cdh.cc <- function(x,par,aargs.par){
  .Call("fpi_mq_cdh",x,par,aargs.par,T,PACKAGE="chifit")
}



