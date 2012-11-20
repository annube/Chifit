

mpi.mq.ob.cc <- function(par,x,aargs){
  .Call("mpi_mq_ob",par,x,aargs,F,PACKAGE="chifit")
}


dmpi.mq.ob.cc <- function(par,x,aargs){
  .Call("mpi_mq_ob",par,x,aargs,T,PACKAGE="chifit")
}

mpi.0.mq.ob.cc <- function(par,x,aargs){
  .Call("mpi_0_mq_ob",par,x,aargs,F,PACKAGE="chifit")
}


dmpi.0.mq.ob.cc <- function(par,x,aargs){
  .Call("mpi_0_mq_ob",par,x,aargs,T,PACKAGE="chifit")
}


fpi.mq.ob.cc <- function(par,x,aargs){
  .Call("fpi_mq_ob",par,x,aargs,F,PACKAGE="chifit")
}


dfpi.mq.ob.cc <- function(par,x,aargs){
  .Call("fpi_mq_ob",par,x,aargs,T,PACKAGE="chifit")
}


