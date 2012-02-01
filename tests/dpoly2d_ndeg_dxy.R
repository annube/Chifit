library(chifit)


par = runif(6)*10 + rnorm(6)
x0 = runif(1,min=-10,max=10)
y0 = runif(1,min=-10,max=10)

##  (1)    1 , x , y , x^2 , xy , y^2

##  d/dx (1) =   0 , 1 , 0 , 2 x , y , 0
dpardx = c( par[2] , 2*par[4] , par[5] , 0 , 0 , 0)

poly2d_ndeg(dpardx,list(xvs=x0,yvs=y0) )

##  d/dy (1) =  0  , 0 , 1 ,  0 , x , 2 y
dpardy = c( par[3] , par[5] , 2 * par[6] , 0, 0, 0)
poly2d_ndeg(dpardy,list(xvs=x0,yvs=y0) )




dpoly2d_ndeg_dxy(par ,list(xvs=x0,yvs=y0) )
