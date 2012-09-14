library(chifit)
library(hadron)

f <- function(x,par) {
  par[1] * log( par[2] * x )
}

df <- function(x,par) {
  cbind(log( par[2] * x ) , par[1]/par[2] )
}


x <- 1:5

par <- c(3.5,7.1 )

y <- f(x,par)
dy <- abs(0.5 + rnorm(length(x))*0.05)
y <- y + dy*rnorm(length(y))


plotwitherror(x,y,dy)
plot(function(x) f(x,par) ,xlim=c(0.5,6),add=T)


so.res <- simple.optim(x,y,dy,f,df,c(-1,1),verbose=T)

wnlls.res <- wnlls(x,y,dy,f.in=f,df.in=df,par=c(1,1))
plot(function(x) f(x,so.res$beta) ,xlim=c(0.5,6),add=T,col="red")


