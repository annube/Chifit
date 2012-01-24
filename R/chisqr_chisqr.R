chisqr_chisqr <- function (pars,ni,breaks,doplot=FALSE){

  nM=sum(ni)
  chisqr=0
  
  theoretical_count = numeric(length(ni))

  n_bins=0
  num_breaks=length(breaks)
  mids = 0.5 * ( breaks[1:(num_breaks-1)]+breaks[2:num_breaks] )

  if(doplot) { 
    print(breaks)
    print(mids)
    print(paste("dof = " , pars[1] ))
  }
  
  for(i in 2:length(breaks)){
   p_e = pchisq( breaks[i]*pars[3],df=pars[1],ncp=pars[2] ) - pchisq(breaks[i-1]*pars[3],df=pars[1],ncp=pars[2] )
   F_i = ni[i-1]
   
   
   
   
   theoretical_count[i] = p_e * nM

  if(doplot) print(theoretical_count[i])

   if( theoretical_count[i] >= 5 #nM*significance 
       ) {
     chisqr=chisqr+(F_i - p_e*nM)^2 / (p_e*nM)
     n_bins=n_bins+1
     if(doplot) points(mids[i-1],theoretical_count[i])
    } else {
      if(doplot) points(mids[i-1],theoretical_count[i],pch=4)
    }
}


  legend("topright",c("theoretical","excluded"),pch=c(as.numeric(1),4))
  return(chisqr)

}



chisqr_test_chisqr <- function (data,chisqr_dof,ncp=0.,stretch=1.0,breaks=20) {

  hist_data <- hist(data,breaks=breaks)

  dof = length(hist_data$counts) - 1

#print( paste ( " df = " , dof ))


  plot(hist_data)

   chisqr_value <- chisqr_chisqr(c(chisqr_dof,ncp,stretch), hist_data$counts ,hist_data$breaks,doplot=TRUE )

  if ( chisqr_value < qchisq(0.95,df=dof) ) { 
# print (paste ( "yeepeeeeeeeeeeeeeeeeeeeeeeeeeee! the data you gave me is very likly gaussian distributed cause:"))
#     print ( paste(" chi^2 = " , chisqr_optim$value , " < qchisqr(0.95,df=dof ) = " , qchisq(0.95,df=dof ) ))
      return (TRUE)
} else   return ( FALSE )


}


chisqr_fit_chisqr <- function (data,chisqr_dof,ncp=0.,stretch=1.0,breaks=20) {

  hist_data <- hist(data,breaks=breaks,plot=FALSE)
  plot(hist_data,ylim=c(min(hist_data$counts),2*max(hist_data$counts)))

  dof = length(hist_data$counts) - 4

#print( paste ( " df = " , dof ))


#plot(hist_data)

   
   optim_res<-optim(par=c(chisqr_dof,ncp,stretch),fn=chisqr_chisqr, ni=hist_data$counts, breaks=hist_data$breaks ,
    method="L-BFGS-B",lower=c(chisqr_dof*0.9,0.0,0.01),upper=c(chisqr_dof*2,30,Inf),control=list(parscale=c(chisqr_dof,1.e-10,1.e-10)))

 chisqr_value <- chisqr_chisqr(optim_res$par, hist_data$counts ,hist_data$breaks,doplot=TRUE )


return(c(optim_res$par,optim_res$value,dof))

}

