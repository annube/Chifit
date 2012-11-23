

mpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_mq_ob_zp",x,par,aargs,F,PACKAGE="chifit")
}


dmpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_mq_ob_zp",x,par,aargs,T,PACKAGE="chifit")
}

mpi.0.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_0_mq_ob_zp",x,par,aargs,F,PACKAGE="chifit")
}


dmpi.0.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("mpi_0_mq_ob_zp",x,par,aargs,T,PACKAGE="chifit")
}


fpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("fpi_mq_ob_zp",x,par,aargs,F,PACKAGE="chifit")
}


dfpi.mq.ob.ZP.cc <- function(x,par,aargs){
  .Call("fpi_mq_ob_zp",x,par,aargs,T,PACKAGE="chifit")
}


