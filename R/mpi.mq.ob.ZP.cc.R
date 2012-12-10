

mpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_mq_ob",x,par,aargs,F,0,T,PACKAGE="chifit")
}


dmpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_mq_ob",x,par,aargs,T,0,T,PACKAGE="chifit")
}

mpi.0.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_0_mq_ob",x,par,aargs,F,0,T,PACKAGE="chifit")
}


dmpi.0.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_0_mq_ob",x,par,aargs,T,0,T,PACKAGE="chifit")
}


fpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("fpi_mq_ob",x,par,aargs,F,0,T,PACKAGE="chifit")
}


dfpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("fpi_mq_ob",x,par,aargs,T,0,T,PACKAGE="chifit")
}


