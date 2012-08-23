

mpi.mq.ob.cc <- function(par,x,aargs){
  .Call("mpi_mq_ob",par,x,aargs,PACKAGE="chifit")
}


dmpi.mq.ob.cc <- function(par,x,aargs){
  .Call("dmpi_mq_ob",par,x,aargs,PACKAGE="chifit")
}
