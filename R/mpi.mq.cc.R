

mpi.mq.cc <- function(par,x,aargs){
  .Call("mpi_mq",par,x,aargs,PACKAGE="chifit")
}


dmpi.mq.cc <- function(par,x,aargs){
  .Call("dmpi_mq",par,x,aargs,PACKAGE="chifit")
}