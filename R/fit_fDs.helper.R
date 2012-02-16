

## general function for fitting and matching of f_Ds


make.poly2d.fit <- function(x1,x2,y,dy,p.degree=2,cache=NULL){
  fit.res <- fit_poly2d(list(xvs=x1,yvs=x2),y,dy,degree=p.degree,confidence=1.0,cache=NULL)
  return( fit.res  )
}


fit.data.poly2d <- function(indices.y,boot.index.y=1,indices.x1=NA,boot.index.x1=NA,indices.x2=NA,boot.index.x2=NA,p.degree=2){
  

  boot.res.y <- extract_boot_results(indices.y,index=boot.index.y)
  boot.R = length(boot.res.y$data_y_boot[,1])

  direct.x1=FALSE
  direct.x2=FALSE
  
  if( ! missing(boot.index.x1) ) {
    direct.x1 = TRUE
    boot.res.x1 <- extract_boot_results(indices.x1,index=boot.index.x1)
    x1 = boot.res.x1$data_y
    x1.boot = boot.res.x1$data_y_boot
  } else {
    x1 = boot.res.y$data_x0
    x1.boot = matrix( rep( x1 , each =  boot.R ) , nrow = boot.R ,ncol = length(x1) )
  }

  if( ! missing(boot.index.x2) ) {
    direct.x2 = TRUE
    boot.res.x2 <- extract_boot_results(indices.x2,index=boot.index.x2)
    x2 = boot.res.x2$data_y
    x2.boot = boot.res.x2$data_y_boot
  } else {
    x2 = boot.res.y$data_x
    x2.boot = matrix( rep( x2 , each =  boot.R ) , nrow = boot.R ,ncol = length(x2) )
  }

  

  fit.res =  make.poly2d.fit(x1,x2,boot.res.y$data_y,boot.res.y$data_dy,p.degree)


  par.boot = array(0,dim=c(boot.R,length(fit.res$par)))
  chisqr.boot = numeric(boot.R)

  if( ! direct.x1 & ! direct.x2 ) {
    cache = fit.res$cache
  } else {
    cache = NULL
  }

  
  for( i in 1:boot.R) {
    
    fit.res.boot =  make.poly2d.fit(x1.boot[i,],x2.boot[i,],boot.res.y$data_y_boot[i,],boot.res.y$data_dy,p.degree,cache=cache)
##    print(fit.res.boot)
    par.boot[i,] = fit.res.boot$par
    chisqr.boot[i] = fit.res.boot$Chisqr

##    print(fit.res.boot$Chisqr)
  }
  

  return(  list(
                fit.res = fit.res,
                chisqr.boot = chisqr.boot,
                par.boot = par.boot,
                p.degree = p.degree
                )
         )

  
  
}

domatch.2d <- function(fit.res,x,x.boot,y,y.boot){
  matching.q <- poly2d_ndeg(fit.res$fit.res$par,list(xvs=x,yvs=y),fit.res$p.degree)

  grad.matching.q <- dpoly2d_ndeg_dxy(fit.res$fit.res$par,list(xvs=x,yvs=y),fit.res$p.degree)
  

  boot.R = length(x.boot)
  
  match.res.boot <- numeric(boot.R)
  for( i in 1:boot.R ){
    match.res.boot[i] = poly2d_ndeg(fit.res$par.boot[i,],list(xvs=x.boot[i],yvs=y.boot[i]),fit.res$p.degree)
  }
  
  return( list(
               matching.q = matching.q,
               grad.matching.q = grad.matching.q,
               match.res.boot = match.res.boot
               )
         )
}


format.result <- function(val,dval,name){
      sprintf( " %s = %f +- %f ",
              name ,
              val ,
              dval
              )
}
