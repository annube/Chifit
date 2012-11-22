library(chifit)



## 
ampi.mq.ob <- function(x,par,aargs){

  R = numeric(aargs$numLs)
  R[1] = 1
  if(aargs$numLs > 1 )
    R[2:aargs$numLs] = par[ 5 + 1:(aargs$numLs-1) ]
  
  print(R)
  
  B <- par[1] * R[aargs$lsIndex]
  f <- par[2] * R[aargs$lsIndex]
  c2 =  par[3] * R[aargs$lsIndex]^4
  Lambda3 <- par[4] * R[aargs$lsIndex]
  CMpm <- par[5] * R[aargs$lsIndex]^2
  M0.sq = 2*B*x + 2 * c2
  Mpm.sq = 2*B*x

  xi.0 <-  M0.sq / ( 4 * pi * f )^2
  
  return(  Mpm.sq * ( 1 +
                      xi.0 * ( log( M0.sq/Lambda3^2) +
                                   g.tilde( sqrt( M0.sq ) * aargs$L , partitions ) ) +
                     CMpm
                     )
         )
  
}

ZPs = c( 0.9 , 0.8 , 0.75 )


x = c(
  0.01,
  0.005,
  0.008,
  0.009,
  0.004,
  0.005)

par = c(
  1.1,
  0.05,
  0.2,
  0.67,
  1.0,
  0.9,
  0.7
       )

lsIndex = c(1,1,2,2,3,3)

aargs = list( numLs = 3,
  lsIndex = lsIndex,
  L = c(24,25,24,48,64,32),
  ZP = ZPs[lsIndex]
  )


mpi.mq.ob.cc(x,par,aargs) -
  ampi.mq.ob(x,par,aargs)


dmpi.mq.ob.cc(x,par,aargs)



par.zp = c(
  1.1,
  0.05,
  0.2,
  0.67,
  1.0,
  0.9,
  0.7,
  ZPs
  )


aargs.zp = list( numLs = 3,
  lsIndex = lsIndex,
  L = c(24,25,24,48,64,32)
  )

par.zp


.Call("mpi_mq_ob_zp",  x , par.zp, aargs.zp, F, PACKAGE = "chifit") - 
mpi.mq.ob.cc(x,par,aargs)



dmpi.mq.ob.cc(x,par,aargs) 
dmpi.mq.ob.ZP.cc(x,par.zp,aargs.zp) 
