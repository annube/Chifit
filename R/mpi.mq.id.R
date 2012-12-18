

mpi.mq.generate.id <- function(quantity,fse) {
  .Call("GenerateChifitExpression",quantity,fse,PACKAGE="chifit")
}

mpi.mq.id <- function(id,x,par,aargs,fitZP=T) {
  .Call("mpi_mq_gen_id", id,
        x,par,aargs,
        F, fitZP,
        PACKAGE="chifit")
}

dmpi.mq.id <- function(id,x,par,aargs,fitZP=T) {
  .Call("mpi_mq_gen_id",id,x,par,aargs,T,fitZP,
        PACKAGE = "chifit" )
}


