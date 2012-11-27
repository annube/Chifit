library(chifit)
mpi.mq.cc(1,c(1,2,3),list(a=1,L=32,ZP=1.) )


dmpi.mq.cc(c(1,2,3,4,5),c(1,2,3),list(a=c(1,2,3,4,5),L=rep(32,5),ZP=rep(1,5)))

fpi.mq.cc(c(1,2,3,4,5),c(1,2,3),list(a=c(1,2,3,4,5),L=rep(32,5),ZP=rep(1,5)))
dfpi.mq.cc(c(1,2,3,4,5),c(1,2,3),list(a=c(1,2,3,4,5),L=rep(32,5),ZP=rep(1,5)))

mpi.mq.ob.cc(c(1,2,3,4,5),c(1,2,3,4),list(a=c(1,2,3,4,5)))
dmpi.mq.ob.cc(c(1,2,3,4,5),c(1,2,3,4),list(a=c(1,2,3,4,5)))
