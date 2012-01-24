


chisqr_norm <- function (pars,ni,breaks,doplot=FALSE){

  nM=sum(ni)
  chisqr=0

  theoretical_count = numeric(length(ni))
  n_bins=0
  num_breaks = length(breaks)
  mids = 0.5 * ( breaks[1:(num_breaks-1)]+breaks[2:num_breaks] )


  for(i in 2:length(breaks)){
   p_e = pnorm( (breaks[i]-pars[1]) / abs( pars[2] ) ) - pnorm( (breaks[i-1]-pars[1]) / abs ( pars[2] ) )
   F_i = ni[i-1]

   theoretical_count[i] = p_e*nM

#   print( paste ( " F_e " ,p_e*nM , " F_i " , F_i , " (F_i - p_e*nM)^2 / (p_e*nM) " ,  (F_i - p_e*nM)^2 / (p_e*nM) ) )
   if( theoretical_count[i] >= 5 #nM*significance 
       ) {
     chisqr=chisqr+(F_i - p_e*nM)^2 / (p_e*nM)
     n_bins=n_bins+1
     if(doplot) points(mids[i-1],theoretical_count[i])
   } else {
     if(doplot) points(mids[i-1],theoretical_count[i],pch='x')
   }
  }

  return( list(chisqr=chisqr,dof = n_bins - 2 ) )

}



chisqr_test_norm <- function (data,breaks=20,doplot=FALSE) {

  mu = mean(data)
  sigma  =  sd(data)
  hist_data <- hist(data,breaks=breaks,plot=doplot)
  

  chisqr.norm <- chisqr_norm(c(mu,sigma) , hist_data$counts ,hist_data$breaks,doplot=doplot )
  chisqr.res <- chisqr.norm$chisqr
  
  dof <- chisqr.norm$dof

  
  if(doplot) print( paste ( " df = " , dof ," chisqr = " , chisqr.res))
  
  if ( chisqr.res < qchisq(0.95,df=dof) ) { 
      # print (paste ( "yeepeeeeeeeeeeeeeeeeeeeeeeeeeee! the data you gave me is very likly gaussian distributed cause:"))
      #     print ( paste(" chi^2 = " , chisqr_optim$value , " < qchisqr(0.95,df=dof ) = " , qchisq(0.95,df=dof ) ))
      return (TRUE)
  } else   return ( FALSE )



}


chisqr_test_optim <- function (data,breaks=20,doplot=FALSE) {

  mu = mean(data)
  sigma  =  sd(data)
  hist_data <- hist(data,breaks=breaks)

  dof = length(hist_data$counts) - 3

  #print( paste ( " df = " , dof ))

  if(doplot)
    plot(hist_data)

  chisqr_optim <- optim( c(mu,sigma) , fn=function(par,x,...) chisqr_norm(par,x,...)$chisqr , ni=hist_data$counts , breaks=hist_data$breaks , method="Nelder-Mead"  )

  xs = seq( min( breaks=hist_data$breaks ) , max( breaks=hist_data$breaks ) , length.out = 100 )

  delta_x = (hist_data$breaks[2] - hist_data$breaks[1])

  ys = sum(hist_data$counts) /sqrt(2 * pi ) / chisqr_optim$par[2] * delta_x * 
             exp( - 0.5 * ( (xs-chisqr_optim$par[1])/chisqr_optim$par[2] )^2 )

  lines ( xs , ys )

  if ( chisqr_optim$value < qchisq(0.95,df=dof) ) { 
      # print (paste ( "yeepeeeeeeeeeeeeeeeeeeeeeeeeeee! the data you gave me is very likly gaussian distributed cause:"))
      #     print ( paste(" chi^2 = " , chisqr_optim$value , " < qchisqr(0.95,df=dof ) = " , qchisq(0.95,df=dof ) ))
      return (TRUE)
  } else   return ( FALSE )


}
