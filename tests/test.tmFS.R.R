library(chifit)

r=0.72
k=2
x0=5


.Call("tmFS_R_fn_R",x0,k,r,PACKAGE="chifit")
.Call("tmFS_R_fn_R_ginac",x0,k,r,PACKAGE="chifit")

pdf("test.tmFS.R.pdf")

## test derivative
f <- function(x) .Call("tmFS_R_fn_R",x,k,r,PACKAGE="chifit")
df <- function(x) .Call("tmFS_R_dx_fn_R",x,k,r,PACKAGE="chifit")


df.approx <- function(x,deri.d=1.e-3) (f(x+deri.d)-f(x))/deri.d


f(x0)
df(x0)

delta = 0.1

plot( function(x) sapply(x, function(y) df(y) ) , xlim=c(x0-delta,x0+delta) )

plot( function(x) sapply(x, function(y) df.approx(y) ) , xlim=c(x0-delta,x0+delta), add=T ,col="red" )

plot( function(x) sapply(x, function(y) df.approx(y,1.e-4) ) , xlim=c(x0-delta,x0+delta), add=T ,col="green" )


plot( function(x) sapply(x, function(y) ( f(y)-f(x0) - df(x0)*(y-x0) ) ) , xlim=c(x0-delta,x0+delta) )



rs = c(0.72,0.82,0.95,1)

plot(  function(x) sapply(
                          x ,
                          function(y) rs[1]^(k+1) * .Call("tmFS_R_fn_R_ginac",y,k,rs[1],PACKAGE="chifit")
                          ),
     xlim=c(2,6)
##     ,log='y'
     )


for( i in 2:4 ) 
plot(  function(x) sapply(
                          x ,
                          function(y) rs[i]^(k+1) * .Call("tmFS_R_fn_R_ginac",y,k,rs[i],PACKAGE="chifit")
                          ),
     xlim=c(2,6),
     add=T,
     col=i
     )

dev.off()
