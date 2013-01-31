
general_fit <- function(which_fit, data_x, data_y, data_dy, y_in , dy_in , x_in , dx_in )
{

  if (which_fit  == "circle" ) {


    fit_function=circle
    dfit_function=dcircle
                                        #  optim_fit <- optim(calc_circle_params(tab1$mu2,tab1$m),
                                        #                	chisqr_mk_function,input_vals=list( fnctn=fit_function,data_x=tab1$mu2,data_y=tab1$m,data_dy=tab1$dm ) ,
                                        #                	method="BFGS")

    print("simulating fit")
    sim_fit <-simulate_fit(parameter = calc_circle_params(data_x,data_y) ,
                           fn = chisqr_mk_function , N=400,
                           input_vals=list( fnctn=fit_function,data_x=data_x,data_y=data_y,data_dy=data_dy ),method="Nelder-Mead" )
    print("simulation done")

    optim_fit <- sim_fit$optim_fit

    if( !missing(y_in) ) {
      ## calculate the matching x value

      par=optim_fit$par

      delta_y = y_in-par[3]
      delta_x = sqrt( par[1] ^ 2 - delta_y ^ 2 )

      x_out_1 = par[2]-delta_x
      x_out_2 = par[2]+delta_x

      
      ## calculate the error

      delta_ys = y_in - sim_fit$simulation[3,]
      
                                        #print(delta_ys)
      
      delta_xs = sqrt( sim_fit$simulation[1,] ^ 2 - delta_ys ^ 2 )
      x_out_1_simulated   = sim_fit$simulation[2,] - delta_xs

      dx_out_1 = sd(x_out_1_simulated)

      ##   
      ##   +
      ## +++++  error from the given y
      ##   +  
      ##

      m = - dcircle( apply (  sim_fit$simulation,1,mean  ) , x_out_1 ) [ 2,1 ]


      dx_out_1 = dx_out_1 + 1/m * ( dy_in )

      return(list(match_y=TRUE,y_in=y_in,dy_in=dy_in,sim_fit=sim_fit,value = x_out_1, dvalue = dx_out_1,bias=mean(x_out_1_simulated)-x_out_1))
    } else if(!missing(x_in)) {
                                        # determine f_K


      
      attach( optim_fit )
      y_out = sqrt( par[1] ^ 2 - (x_in - par[2] )^2 ) + par[3]
      detach( optim_fit )
      
                                        # determine error of f_K
      
      
      attach( sim_fit)
      y_outs = sqrt( simulation[1,] ^ 2 - (x_in - simulation[2,] )^2 ) + simulation[3,]
      
      dy_out = sd(y_outs)
      
      m_y_out = - dcircle( apply (  simulation,1,mean  ) , x_in ) [ 2,1 ]
      
      detach( sim_fit ) 
      
                                        # ++ error from m_s
      
      dy_out = dy_out + m_y_out * dx_in

      return(list(match_y=FALSE,x_in=x_in,dx_in=dx_in,sim_fit=sim_fit,value = y_out, dvalue = dy_out,bias=mean(y_outs)-y_out))
    } else {return (sim_fit) }
    

  }
}

general_fit_boot <- function(which_fit, boot_results, y_in , dy_in , x_in , dx_in,poly_degree, ...){

  ## contains data_x, data_y, data_dy, data_y_boot
  attach(boot_results) 

  ## determine the number of boot samples
  num_boot_samples = length(data_y_boot[,1])


  ## so far only circle is implemented


  print("simulating fit via boot")


  ## find NaN's
  sel_range_l = ! is.nan(data_y_boot[,1])

  ## find the largest subset among all data sets (belonging to different masses) containing no Nan
  for(j in 2:length(data_y_boot[1,]))
    sel_range_l = sel_range_l & ! is.nan(data_y_boot[,j])


  ##detect outliers
  
  ## eventually not Nan have been detected directly, it was observed that Nan in earlier stages of
  ## calculations lead to outliers in later stages
  ## I have implemented a little outlier-detection algorithm to sort out these outliers
  for(j in 1:length(data_y_boot[1,])) {
    outly <- outliers(data_y_boot[,j])

                                        # find again the largest subset containing no Nan's and no outliers
    sel_range_l = sel_range_l & ! outly 
  }

                                        # convert again to indices
  sel_range =  which(sel_range_l)
  



  print(which(!sel_range_l))


  if (which_fit  == "circle" ) {

    sim_optim_fit = array(NA,dim=c(num_boot_samples,4))
    
    ## set fit functions
    fit_function=circle
    dfit_function=dcircle

    
    for( i in 1:num_boot_samples ) {

      if(sel_range_l[i]){
        ##         print(i)
        par = calc_circle_params(data_x,data_y_boot[i,])
        optim_res <- optim(par=par,
                                   fn=chisqr_mk_function,gr=dchisqr_mk_function,
                                   input_vals=list( fnctn=fit_function,dfnctn=dfit_function,data_x=data_x,data_y=data_y_boot[i,],data_dy=data_dy ),method="Nelder-Mead")
        sim_optim_fit[i,1:3] <- optim_res$par
        sim_optim_fit[i,4] <- optim_res$value
        
      }
    }
    
    
    print("simulation done")
    
    optim_fit <- optim(par=calc_circle_params(data_x,data_y),
                       fn=chisqr_mk_function,
                       input_vals=list( fnctn=fit_function,data_x=data_x,data_y=data_y,data_dy=data_dy ),method="Nelder-Mead")
    
  } else if ( which_fit == "poly" ) {
    fit_function=mypoly
    dfit_function=dmypoly
    

    if( missing(poly_degree) )
      poly_degree <- 2
    
    npars <- (poly_degree+1)

    sim_optim_fit = array(NA,dim=c(num_boot_samples,npars+2))

    
    for( i in 1:num_boot_samples ) {

      if(sel_range_l[i]){
       
        fit_tmp<-fit_multilinear_function(data_x,data_y_boot[i,],
                                          data_dy,
                                          fit_function,
                                          dfit_function,
                                          npars=npars,
                                          maxdeg=poly_degree,
                                          ...)
        sim_optim_fit[i,1:npars] <- fit_tmp$par
        sim_optim_fit[i,(npars+1)] <- fit_tmp$Chisqr
        sim_optim_fit[i,(npars+2)] <- fit_tmp$unused
      }
    }
    
    optim_fit <-fit_multilinear_function(data_x,data_y,
                                                   data_dy,
                                                   fit_function,
                                                   dfit_function,
                                                   npars=npars,
                                                   maxdeg=poly_degree)

      
    
  }


  
  
  detach(boot_results)

  input_vals=boot_results
  input_vals[["fnctn"]] <- fit_function



  sim_fit = list(input_vals=input_vals,optim_fit=optim_fit,sim_optim_fit=sim_optim_fit)


  ## go to x matching mode
  if( !missing(y_in) ) {
    ## calculate the matching x value

    ## take parameter from fit
    par=optim_fit$par

    if( which_fit == "circle" ) {
      delta_y = y_in-par[3]
      delta_x = sqrt( par[1] ^ 2 - delta_y ^ 2 )

      ## calculate matching x
      x_out_1 = par[2]-delta_x
      x_out_2 = par[2]+delta_x

      x_out <- x_out_1

    
      ## calculate the error via (manual) boot

      delta_ys = y_in - sim_optim_fit[,3]
      delta_xs = sqrt( sim_optim_fit[,1] ^ 2 - delta_ys ^ 2 )
      
      x_out_1_simulated   = sim_optim_fit[,2] - delta_xs
      x_out_2_simulated   = sim_optim_fit[,2] + delta_xs


      x_out_simulated = x_out_1_simulated
      
      print(x_out_simulated[sel_range])
      print(x_out_simulated[! sel_range_l])
    } else if( which_fit == "poly" ) {
      print("not implemented so far, and there is so far no reason for doing it
as the circular fit meets all requirements" )
    }


    ## create boot object for later use with boot.ci
    boot_x_out = list()
    boot_x_out$t=array(x_out_simulated[sel_range],dim=c(length(sel_range),1))
    boot_x_out$R=length(sel_range)
    boot_x_out$t0=x_out
    boot_x_out$call="boot"
    attr(boot_x_out,"class") <- "boot"

    x_out_boot.ci <- boot.ci(boot_x_out,type="norm")

    ## rather use boot estimate of confidence interval 
    dx_out = (x_out_boot.ci$normal[3]-x_out_boot.ci$normal[2])/2/1.96

    ##	dx_out = sd(x_out_1_simulated)
    print(paste("error from boot ", dx_out) )
    ##   
    ##   +
    ## +++++  error from the given y
    ##   +  
    ##

    m = abs(dcircledx( optim_fit$par , x_out ))

    print(paste(" slope of fitting circle in fitting region is " , m ))


    dx_out = dx_out + 1/m * ( dy_in )

    dx_out_ppe <- 1/m * ( dy_in )

    print(paste("error from propagated error ", 1/m * ( dy_in ) ) )


    return(list(match_y=TRUE,y_in=y_in,dy_in=dy_in,sim_fit=sim_fit,value = x_out, dvalue = dx_out, dvalue_ppe = dx_out_ppe,
                x_out_simulated=x_out_simulated,sel_range_l=sel_range_l,
                bias=mean(x_out_simulated[sel_range])-x_out))
  } else if(!missing(x_in)) {
                                        # determine f_K

    par=optim_fit$par

    y_out = sqrt( par[1] ^ 2 - (x_in - par[2] )^2 ) + par[3]
    print(paste("par[] " , par[1] , par[2] ,  par[1] , " x_in" , x_in , "y_out" , y_out ) )
    
    ## determine error of f_K
    
    
    y_outs = apply(sim_optim_fit,1,circle,x_in)

    print(y_outs)
    ##print(sim_optim_fit[sel_range_l,1])
    ##print(sim_optim_fit[! sel_range_l,1])

    
                                        # create boot object for later use with boot.ci
    boot_y_out = list()
    boot_y_out$t=array(y_outs[sel_range],dim=c(length(sel_range),1) )
    boot_y_out$R=length(sel_range)
    boot_y_out$t0=y_out
    boot_y_out$call="boot"
    attr(boot_y_out,"class") <- "boot"
    y_out_boot.ci <- boot.ci(boot_y_out,type="norm")

                                        # rather use boot estimate of confidence interval 
    dy_out = (y_out_boot.ci$normal[3]-y_out_boot.ci$normal[2])/2/1.96

                                        #dy_out = sd(y_outs)
    
    print(paste("error from boot ",dy_out) )
    
    m_y_out = abs( dcircledx( optim_fit$par , x_in ) )

                                        # ++ error from m_s
    dy_out = dy_out + m_y_out * dx_in

    print(paste("error from propagation ", m_y_out * dx_in) )


    return(list(match_y=FALSE,x_in=x_in,dx_in=dx_in,sim_fit=sim_fit,value = y_out, dvalue = dy_out,
                y_out_simulated=y_outs,sel_range_l=sel_range_l,
                bias=mean(y_outs)-y_out))
  } else {return (sim_fit) }
  


}



direct_fit_x <- function(boot_x,boot_y,y_in,dy_in){

  num_boot_samples = length(boot_x$data_y_boot[,1])

  lm_fit <- lm(boot_y$data_y~boot_x$data_y)


  m=coefficients(lm_fit)[2]
  n=coefficients(lm_fit)[1]

  x_out = (y_in-n)/m
  dx_out = dy_in/m


  print("simulating fit via boot")
  

  
  lm_fit_boot = array(dim=c(num_boot_samples,3))

  for( i in 1:num_boot_samples) {
    lm_fit_tmp <- lm(boot_y$data_y_boot[i,]~boot_x$data_y_boot[i,],weights=1/boot_y$data_dy^2)
    lm_fit_boot[i,1:2] <- lm_fit_tmp$coefficients
    lm_fit_boot[i,3]=sum(lm_fit_tmp$residuals^2)
  }
  
  print("simulation done")


  x_out_boot = (y_in - lm_fit_boot[,1] )/ lm_fit_boot[,2]
  dx_out = dx_out + sd(x_out_boot)

  return(list(data_x=boot_x$data_y,data_dx=boot_x$data_dy,
              data_y=boot_y$data_y,data_dy=boot_y$data_dy,
              y_value=y_in,dy_value=dy_in,
              x_value=x_out,dx_value=dx_out,
              lm_fit = lm_fit,
              lm_fit_boot=lm_fit_boot,
              boot_results_x=boot_x,
              x_out_boot = x_out_boot, y_out_boot = rep(y_in,num_boot_samples),
              boot_results_y = boot_y, swap_xy=FALSE))

}


direct_fit_y <- function(boot_x,boot_y,x_in,dx_in){

  num_boot_samples = length(boot_x$data_y_boot[,1])

  lm_fit <- lm(boot_y$data_y~boot_x$data_y)


  m=coefficients(lm_fit)[2]
  n=coefficients(lm_fit)[1]

  y_out = m*x_in + n
  dy_out = dx_in * m


  print("simulating fit via boot")
  

  
  lm_fit_boot = array(dim=c(num_boot_samples,3))

  for( i in 1:num_boot_samples) {
    lm_fit_tmp <- lm(boot_y$data_y_boot[i,]~boot_x$data_y_boot[i,],weights=1/boot_y$data_dy^2)
    lm_fit_boot[i,1:2] <- lm_fit_tmp$coefficients
    lm_fit_boot[i,3]=sum(lm_fit_tmp$residuals^2)
  }
  
  print("simulation done")


  y_out_boot = x_in * lm_fit_boot[,2] + lm_fit_boot[,1]
  dy_out = dy_out + sd(y_out_boot)

  return(list(data_x=boot_x$data_y,data_dx=boot_x$data_dy,
              data_y=boot_y$data_y,data_dy=boot_y$data_dy,
              y_value=y_out,dy_value=dy_out,
              x_value=x_in,dx_value=dx_in,
              lm_fit=lm_fit,
              lm_fit_boot=lm_fit_boot,
              boot_results_x=boot_x,
              x_out_boot = rep(x_in,num_boot_samples) , y_out_boot = y_out_boot,
              boot_results_y = boot_y,swap_xy = TRUE ))

}



                                        # $Id: plotutils.R 134 2010-05-13 16:06:21Z urbach $
plotwitherror_dx <- function(x, dx=NULL , y, dy,xlim ,ylim, rep=FALSE, ... ,arrow.l=0.02) {

  if(missing(ylim))   ylim=c(min(y-2*dy, na.rm = TRUE),max(y+2*dy, na.rm = TRUE))
  if(missing(xlim))   xlim=c(min(x-2*dx, na.rm = TRUE),max(x+2*dx, na.rm = TRUE))

  if(!rep){
    plot(x,y, xlim = xlim, ylim = ylim,...)
  } else {
    points(x,y, xlim = xlim, ylim = ylim,...)
  }
    
  arrows(x, y-dy, x, y+dy, length=arrow.l,angle=90,code=3,... )
  if( ! missing( dx ) )
    arrows(x-dx, y, x+dx, y, length=arrow.l,angle=90,code=3,... )
}

plot_direct_fit <- function(fit_result,ensemble_name="YZZ.VV",meson="X"){

  pdf(file=paste(ensemble_name,"_plot_m_vs_f_",meson, ".pdf",sep=""))
  attach(fit_result)

  x_min=min(c(data_x-2*data_dx,x_value-2*dx_value))
  x_max=max(c(data_x+2*data_dx,x_value+2*dx_value))

  y_min=min(data_y-2*data_dy)
  y_max=max(data_y+2*data_dy)

  labels_pre = c(bquote(f[.(meson)]),bquote(m[.(meson)]))

  if(swap_xy)
    {labels = labels_pre[c(2,1)]}
  else 
    { labels = labels_pre }

  plotwitherror_dx(data_x,data_dx,data_y,data_dy,xlab=labels[1],ylab=labels[2],ylim=c(y_min,y_max), xlim=c(x_min,x_max) )


  lines( c(x_min,x_value) , rep( y_value-dy_value,2 ),lty="dashed" )
  lines( c(x_min,x_value) , rep( y_value,2 ) )
  lines( c(x_min,x_value) , rep( y_value+dy_value,2 ),lty="dashed" )

  lines( rep(x_value-dx_value,2) , c(y_min,y_value) ,lty="dashed" )
  lines( rep(x_value,2) , c(y_min,y_value)  )
  lines( rep(x_value+dx_value,2) , c(y_min,y_value) ,lty="dashed" )



  xs = c(x_min,x_max)

  lines(xs, lm_fit$coefficients[2]*xs + lm_fit$coefficients[1] )

  for( i in 1:10 ) {
    lines(xs, lm_fit_boot[i,2]*xs + lm_fit_boot[i,1] , col="gray")
    points (boot_results_x$data_y_boot[i,],boot_results_y$data_y_boot[i,],pch=3,col="gray")
  } 

  points(x_out_boot,y_out_boot,pch=4)

  text(  0.5*(x_min+x_max), y_min ,  bquote( chi^2 == .( sum(lm_fit$residuals^2)  ) ) ) 


  hist(lm_fit_boot[,3])

  detach(fit_result)
  dev.off()

}


plot_fit <- function(sim_fit_,ensemble_name="YZZ.VV",meson="X",force_plot = FALSE) {
  attach(sim_fit_)
  attach(sim_fit_$sim_fit$input_vals)

  if( meson == "K" ) {
    quark = "s"
  } else if (meson == "D" ) {
    quark = "c"
  } else if( meson == "X" ) {
    quark = "z"
  }


  
  x_min=min(data_x)
  x_max=max(data_x)
  
  y_min=min(data_y)
  y_max=max(data_y)
  
  y_min_f=min(data_y-2*data_dy)
  y_max_f=max(data_y+2*data_dy)

  if( match_y ) {

    x_min=min(x_min,value-2*dvalue)
    x_max=max(x_max,value+2*dvalue)
    
    y_min_f=min(y_min_f,y_in-2*dy_in)
    y_max_f=max(y_max_f,y_in+2*dy_in)


  }
  

  xs = seq(x_min,x_max,length.out=500)


  fit_fn_vals = fnctn(sim_fit_$sim_fit$optim_fit$par,xs)


  if( match_y ) {

    pdffile = paste(ensemble_name,"_plot_m_",meson, ".pdf",sep="")
    if( file.exists(pdffile) & ! force_plot) {
        detach(sim_fit_$sim_fit$input_vals)
        detach(sim_fit_)
        return()
    }

    pdf(file=pdffile)
    

    plotwitherror(data_x,data_y,data_dy,pch=12,xlab=bquote( m[.(quark)] ) , ylab=bquote(m[.(meson)]),xlim=c(x_min,x_max),ylim=c(y_min_f,y_max_f))

    lines(xs,fit_fn_vals)





    ## plot boot envelope 
    ys_fit <- apply( sim_fit$sim_optim_fit[,],1,fnctn,value)

    points(array(value,dim=c(length(ys_fit))),ys_fit,pch="-")



    max_ys_fit=fnctn(sim_fit$sim_optim_fit[1,],xs)
    min_ys_fit=fnctn(sim_fit$sim_optim_fit[1,],xs)

    for(i in 1: length(ys_fit)){
      max_ys_fit = apply(data.frame(fnctn(sim_fit$sim_optim_fit[i,],xs),max_ys_fit),1,max)
      min_ys_fit = apply(data.frame(fnctn(sim_fit$sim_optim_fit[i,],xs),min_ys_fit),1,min)
    }

    lines(xs,max_ys_fit,lty="dashed",col="grey")
    lines(xs,min_ys_fit,lty="dashed",col="grey")



    alphax = 0.3
    posx = (1-alphax)*x_min + alphax*value
    text(  posx , y_in+dy_in  ,   bquote(paste(m[.(meson)]^{unitary} == .(y_in)," (Input) " ))   ,pos=3 ) 
    
    lines( c(x_min,value) , rep( y_in-dy_in,2 ),lty="dashed" )
    lines( c(x_min,value) , rep( y_in,2 ) )
    lines( c(x_min,value) , rep( y_in+dy_in,2 ),lty="dashed" )
    
    
    
    if((value - x_min)/(x_max-x_min) < 0.5){
      positioning = 4
    } else{
      positioning = 2
    }
    
    alphay = 0.5
    posy = (1-alphay)*y_in+ alphay*y_min
    text(  value+dvalue, posy ,  bquote(m[.(quark)] == .(value) ) ,pos=positioning) 

    ## plot the chi square value
    text(  0.5*(x_min+x_max), y_min-2*data_dy[1] ,  bquote( chi^2 == .(sim_fit$optim_fit$value) ) ) 
    
    
    lines( rep(value-dvalue,2) , c(y_min,y_in) ,lty="dashed" )
    lines( rep(value,2) , c(y_min,y_in)  )
    lines( rep(value+dvalue,2) , c(y_min,y_in) ,lty="dashed" )
    
    title(main=bquote(paste(m[.(quark)] ," from " , m[.(meson)]^{unitary} == m[.(meson)]^{OS} , " Ensemble: " , .(ensemble_name))))

    
    dev.off()
  } else {

    pdffile = paste( ensemble_name,"_plot_f_",meson,".pdf",sep="")
    if(file.exists(pdffile) & ! force_plot) {
      detach(sim_fit_$sim_fit$input_vals)
      detach(sim_fit_)
      return()
    }
    pdf(file=pdffile)
    
    plotwitherror(data_x,data_y,data_dy,pch=12,xlab=bquote( m[.(quark)] ) , ylab=bquote(f[.(meson)]) )
    

                                        # plot boot envelope 
    ys_fit <- apply( sim_fit$sim_optim_fit[,],1,fnctn,value)

    points(array(value,dim=c(length(ys_fit))),ys_fit,pch="-")



    max_ys_fit=fnctn(sim_fit$sim_optim_fit[1,],xs)
    min_ys_fit=fnctn(sim_fit$sim_optim_fit[1,],xs)

    for(i in 1: length(ys_fit)){
      max_ys_fit = apply(data.frame(fnctn(sim_fit$sim_optim_fit[i,],xs),max_ys_fit),1,max)
      min_ys_fit = apply(data.frame(fnctn(sim_fit$sim_optim_fit[i,],xs),min_ys_fit),1,min)
    }

    lines(xs,max_ys_fit,lty="dashed",col="grey")
    lines(xs,min_ys_fit,lty="dashed",col="grey")


    lines(xs,fit_fn_vals)
    
    alphax = 0.3
    posx = (1-alphax)*x_min + alphax*x_in 
    text(  posx , value+dvalue  ,   bquote(f[.(meson)] == .(value) )   ,pos=3 ) 
    
    lines( c(x_min , x_in) , rep( value+dvalue,2),lty="dashed" )
    lines( c(x_min , x_in) , rep( value,2) )
    lines( c(x_min , x_in) , rep( value-dvalue,2),lty="dashed" )
    
                                        #arrows( m_s_1 , f_K, x_min,f_K)
    
    
    if((x_in - x_min)/(x_max-x_min) < 0.5){
      positioning = 4
    } else{
      positioning = 2
    }
    
    alphay = 0.5
    posy = (1-alphay)*value + alphay*y_min_f
    text(  x_in+dx_in, posy ,  bquote(paste(m[.(quark)] == .(x_in) , " (Input) " )) ,pos=positioning) 

    ## plot the chi square value
    text(  0.5*(x_min+x_max), y_min-2*data_dy[1] ,  bquote( chi^2 == .(sim_fit$optim_fit$value) ) ) 
    
    lines( rep(x_in-dx_in,2) , c(y_min_f,value) ,col="darkgreen" ,lty="dashed" )
    lines( rep(x_in,2) , c(y_min_f,value) ,col="darkgreen" )
    lines( rep(x_in+dx_in,2) , c(y_min_f,value) ,col="darkgreen",lty="dashed" )
                                        #axis(1,at=m_s_1,cex.axis=0.8)
    
    
                                        #arrows( m_s_1 , f_K , m_s_1, y_min_f ,col="darkgreen" )
    
    title(main=bquote(paste(f[.(meson)] ," from " ,  {f[.(meson)]^{OS}} ( m[.(quark)] ) , " Ensemble: " , .(ensemble_name) )))
    
    
    dev.off()

  }

  detach(sim_fit_$sim_fit$input_vals)
  detach(sim_fit_)



}
