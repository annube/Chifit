


residue.vector <- function(par,x) {
  d = x - par[3:4]
  v= par[1:2]
  return ( v * ( v %*% d ) / sum(v^2) - d) 
}


distance.chisqr <- function(par,x,dx,y,dy) {
  Rs<-apply(rbind(x,y),2,function(x) residue.vector(par,x) )
#   points(x+Rs[1,],y+Rs[2,])
#   arrows(x, y, x+Rs[1,], y+Rs[2,], length=0.01,angle=90,code=3 )
  chisqr=sum( (Rs[1,]/dx)^2 ) + sum( (Rs[2,]/dy)^2 )
  return(chisqr)
}


distance.regression <- function(x,dx,y,dy) {

  par=c(1,1,0,0)
  optim.res <- optim(par,distance.chisqr,x=x,dx=dx,y=y,dy=dy,method="BFGS")

  par = optim.res$par

  m=par[2]/par[1]

  n=-par[2]/par[1]*par[3]+par[4]

  return(list(m=m,n=n,chisqr=optim.res$value))

}