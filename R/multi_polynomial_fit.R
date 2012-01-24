

mypoly <- function(par,data_x,maxdeg){

  if(missing(maxdeg))
    maxdeg=length(par)-1
  
  ys = rep(par[maxdeg+1],length(data_x))
 
  for( i in maxdeg:1 ) {
    ys = ys*data_x + par[i]
  }

  return( ys )
}

dmypoly <- function(par,data_x,maxdeg = 2){

  result <- array(0,dim=c((maxdeg+1),length(data_x) ) )

  data_x_n = rep(1,length(data_x) )

  for( i in 0:maxdeg ) {
    result[(i+1),] = data_x_n
    data_x_n = data_x_n * data_x
  }

  return( result )
  
}

dmypolydx <- function(par,data_x,maxdeg = 2){
  
  ys = rep(par[maxdeg+1]*maxdeg,length(data_x))
 
  for( i in maxdeg:2 ) {
    ys = ys*data_x + par[i]*(i-1)
  }

  return( ys )
}

dmypolydx_par<- function(par){

  par_ret=par[2:length(par)]
  par_ret=par_ret*1:length(par_ret)

  return(par_ret)
  
}




mypoly_as_fn <- function(par,maxdeg=2) {
  return( function(x) mypoly(par,x,maxdeg) )
}

mypoly_from_roots <- function(roots,data_x,C=1) {
  ys = rep(C,length(data_x))
  for( i in 1:length(roots) )
    ys = ys * ( data_x - roots[i] )

  return(ys)
}


poly2d_ndeg <- function(par,data_x,maxdeg=2){

  attach(data_x)
  
  n_points=length(xvs)

  result <- numeric(n_points)
  count = 1
  for( i in 0:maxdeg) {  
    for( j in 0:i)  {
      result=result + par[count]*xvs^(i-j)*yvs^(j)
      count = count + 1
    }
  }
  
  detach(data_x)

  return(result)
}


dpoly2d_ndeg <- function(par,data_x,maxdeg=2){

  npars = ((maxdeg+1)*(maxdeg+2))/2

  attach(data_x)

  result <- array(0,dim = c(npars,length(xvs)))
  count = 1
  for( i in 0:maxdeg) {  
    for( j in 0:i)  {
      result[count,]=xvs^(i-j) * yvs^(j)
      count = count + 1
    }
  }
  
  detach(data_x)

  return(result)
}



fit_poly2d <-function( data_x,data_y,data_dy,degree=2,cache=NULL,...){

  npars=((degree+1)*(degree+2))/2
  res <-          fit_multilinear_function(
                                  data_x,data_y,
                                  data_dy,
                                  fn=poly2d_ndeg,dfn=dpoly2d_ndeg,
                                  npars=npars,
                                  maxdeg=degree,cache=cache,...)

  res[["degree"]] <- degree
  
  return( res  )
  
}

fit_poly <-function( data_x,data_y,data_dy,degree=2,cache=NULL,...){

  npars=(degree+1)
  res <-          fit_multilinear_function(
                                  data_x,data_y,
                                  data_dy,
                                  fn=mypoly,dfn=dmypoly,
                                  npars=npars,
                                  maxdeg=degree,cache=cache,...)

  res[["degree"]] <- degree
  
  return( res  )
  
}

fit_multilinear_function<-function( data_x,data_y,data_dy,fn,dfn,npars,...,cache=NULL,par,confidence=0.1){

  if( missing(par) ) {
    par=numeric(npars)
  }
  

  if( missing(data_dy)  )
    data_dy = rep(1,length(as.vector(data_y)))
  
  if( is.null(cache) ) {
  
  cache=list()
  
  # retrieve poly (non orthogonal) basis vector
  cache$V = dfn(par,data_x,...)
#  print(dim(cache$V))
  

#  print("ortho")
  # orthogonalize them
  qr_V <- qr(t(cache$V))
  cache$W <- qr.Q(qr_V)
  
#  print("get trafo")
# get transformation matrix
  cache$R <- qr.R(qr_V)
  
  }
  
  attach(cache)

#  print(dim(W))
  
  # perform linear regression
  lm_fit <- lm(as.vector(data_y)~W+0,weights=1/data_dy^2)
  sum_lm_fit<-summary(lm_fit)
  
  # get coefficients in orthonormal basis
  par<-sum_lm_fit$coefficients[,1]
  coeff_dont_use <- which( sum_lm_fit$coefficients[,4] > confidence )
  
  # neglect non-significant ones
 par[coeff_dont_use]=0

##  if (sum( sum_lm_fit$coefficients[,4] > confidence  )>0) print(coeff_dont_use )
  
  # transform back to initial basis
  C = solve(R,par)


  z_hat = fn(C,data_x,...)
 
  Chisqr = sum( ( (z_hat-as.vector(data_y))/data_dy ) ^2)

  detach(cache)

  return(list(par=C,data_x=data_x,data_y=data_y,data_dy=data_dy,unused=sum(coeff_dont_use),Chisqr=Chisqr,cache=cache,lm_fit=lm_fit) )

}


make_grid <- function(x,y){
nx = length(x)
ny = length(y)

xvs = numeric(nx*ny)
yvs = numeric(nx*ny)

attr(xvs,"dim") <- c(nx,ny)
attr(yvs,"dim") <- c(nx,ny)

for( i in 1:ny )
  xvs[,i] = x
for( i in 1:nx )
  yvs[i,] = y

attr(xvs,"dim") <- c(nx*ny)
attr(yvs,"dim") <- c(nx*ny)
return(list(xvs=xvs,yvs=yvs))
}


plot_poly_fit <- function(poly_fit,file="rgl.plot.png",surface_res = 10,scale=c(2,1,0.25),radius=0.0005,write_file=FALSE,theta=-120,phi=15){
  library(rgl)

  
  xs <- seq(min(poly_fit$data_x$xvs),max(poly_fit$data_x$xvs),length.out=surface_res)
  ys <- seq(min(poly_fit$data_x$yvs),max(poly_fit$data_x$yvs),length.out=surface_res)
  
  grid_xy <- make_grid(xs,ys)
  z_hat = poly2d_ndeg(poly_fit$par,poly_fit$data_x,maxdeg=poly_fit$degree)
  z_surface = poly2d_ndeg(poly_fit$par,grid_xy,maxdeg=poly_fit$degree)
  
  attr(z_surface,"dim") <- c(surface_res, surface_res)
  
  rgl.open()

   
  rgl.viewpoint(scale=scale)
  rgl.spheres(poly_fit$data_x$xvs,
              poly_fit$data_y,
              #z_hat,
              poly_fit$data_x$yvs,radius=radius,color="green")
  if(! is.null(poly_fit$fit_point) ){
    rgl.spheres(poly_fit$fit_point[1],poly_fit$fit_point[3],poly_fit$fit_point[2],radius=radius,color="green")
  }
  
  
  rgl.surface( xs, ys, z_surface,alpha=c(0.3),color="red")
  rgl.surface( xs, ys, z_surface,front="line",back="line",alpha=c(0.5),color="black")
  

  
  rgl.bbox()
  
  title3d("Fitting Result")
  
  rgl.viewpoint(theta=theta,phi=phi)
  
  if(write_file) rgl.postscript(filename=paste(file,".pdf",sep=""),fmt="pdf")
  if(write_file) rgl.snapshot(filename=paste(file,".png",sep=""),fmt="png")
  
  rgl.open()
  rgl.viewpoint(theta=theta,phi=phi,scale=scale*c(1,100,1))

  z_surface[,]=0  
  rgl.surface( xs, ys, z_surface,alpha=c(0.3),color="red")
  rgl.surface( xs, ys, z_surface,front="line",back="line",alpha=c(0.5),color="black")
  rgl.spheres(poly_fit$data_x$xvs,poly_fit$data_y-z_hat,poly_fit$data_x$yvs,radius=radius,color="green")
  rgl.bbox()
  
  title3d("Residues")
# if(write_file) rgl.postscript(file,fmt="pdf")


}



test_poly_fit <- function() {

N = 10


degree = 4
npars = ((degree+1)*(degree+2))/2

x = seq(-1,1,length.out=N)
y = seq(-1,1,length.out=N)

nx = length(x)
ny = length(y)

xy_grid <- make_grid(x,y)

xvs=xy_grid$xvs
yvs=xy_grid$yvs



z = sin( 4*pi^2*x%o%y )

R=10
x0=0
y0=0

for( i in 1:nx )
for( j in 1:ny)
 z[i,j] = sqrt(R^2-(x[i]-x0)^2-(y[j]-y0)^2)

dz = array(1,dim=c(nx,ny))

#fit_poly2d <-function( data_x,data_y,data_dy,degree=2){

C = fit_poly2d(data_x=list(xvs=xvs,yvs=yvs),
           data_y = as.vector(z),
           data_dy = rep(1,length(as.vector(z))),degree=degree)


# apply polynomial
z_hat = poly2d_ndeg(C,list(xvs=xvs,yvs=yvs),maxdeg=degree)

print(sum(abs(z_hat-as.vector(z))))

}

#test_poly_fit()
