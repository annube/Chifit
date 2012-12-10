

mpi.mq.ob.cc <- function(x,par,aargs,FSE=0,fitZP=F){
  .Call("mpi_mq_ob",x,par,aargs,F,FSE,fitZP,PACKAGE="chifit")
}


dmpi.mq.ob.cc <- function(x,par,aargs,FSE=0,fitZP=F){
  .Call("mpi_mq_ob",x,par,aargs,T,FSE,fitZP,PACKAGE="chifit")
}

mpi.0.mq.ob.cc <- function(x,par,aargs,FSE=0,fitZP=F){
  .Call("mpi_0_mq_ob",x,par,aargs,F,FSE,fitZP,PACKAGE="chifit")
}


dmpi.0.mq.ob.cc <- function(x,par,aargs,FSE=0,fitZP=F){
  .Call("mpi_0_mq_ob",x,par,aargs,T,FSE,fitZP,PACKAGE="chifit")
}


fpi.mq.ob.cc <- function(x,par,aargs,FSE=0,fitZP=F){
  .Call("fpi_mq_ob",x,par,aargs,F,FSE,fitZP,PACKAGE="chifit")
}


dfpi.mq.ob.cc <- function(x,par,aargs,FSE=0,fitZP=F){
  .Call("fpi_mq_ob",x,par,aargs,T,FSE,fitZP,PACKAGE="chifit")
}


