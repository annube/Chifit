

mpi.mq.lso.cc <- function(x,par,aargs.par){
  .Call("mpi_mq_lso",x,par.par,aargs.par,F,PACKAGE="chifit")
}


dmpi.mq.lso.cc <- function(x,par,aargs.par){
  .Call("mpi_mq_lso",x,par.par,aargs.par,T,PACKAGE="chifit")
}

fpi.mq.lso.cc <- function(x,par,aargs.par){
  .Call("fpi_mq_lso",x,par,aargs.par,F,PACKAGE="chifit")
}


dfpi.mq.lso.cc <- function(x,par,aargs.par){
  .Call("fpi_mq_lso",x,par,aargs.par,T,PACKAGE="chifit")
}
