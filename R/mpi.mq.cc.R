

mpi.mq.cc <- function(par,x,aargs){
  .Call("mpi_mq",par,x,aargs,PACKAGE="chifit")
}


dmpi.mq.cc <- function(par,x,aargs){
  .Call("dmpi_mq",par,x,aargs,PACKAGE="chifit")
}

fpi.mq.cc <- function(par,x,aargs){
  .Call("fpi_mq",par,x,aargs,PACKAGE="chifit")
}


dfpi.mq.cc <- function(par,x,aargs){
  .Call("dfpi_mq",par,x,aargs,PACKAGE="chifit")
}

