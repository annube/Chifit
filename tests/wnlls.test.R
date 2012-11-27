library(chifit)
expfn<-function(x,par,...) par[1] * exp(-par[2]*x)
dexpfn<-function(x,par,...) cbind(  exp(-par[2]*x) , -par[1]*x*exp(-par[2]*x) )
x=0:32
par=c(0.1234,.2)
ys=expfn(x,par)
dys=ys*0.1
plotwitherror(x,ys,dys)
plotwitherror(x,ys,dys,log='y')
ys=expfn(x,par)*( 1 + 0.1 /2 * rnorm(33) )
plotwitherror(x,ys,dys,log='y')
C = diag(dys^2)


wnlls.res <- wnlls(x.in=x,y=ys,
##                   dy=dys,
                   C=C,
                   f=expfn     ,
                   df = dexpfn,
                   par=c(0.1,0.5),
##                   range=10:20
                   )


plot.wnlls(wnlls.res)
