library(chifit)

r=0.72
k=2
x0=5


.Call("tmFS_R_fn_R",x0,k,r,PACKAGE="chifit")
.Call("tmFS_R_fn_R_ginac",x0,k,r,PACKAGE="chifit")

pdf("test.tmFS.R.pdf")


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
