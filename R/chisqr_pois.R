


chisqr_pois <- function (pars,ni,mids,significance=0.05){

  nM=sum(ni)
  chisqr=0

  theoretical_count = numeric(length(ni))
  n_bins=0

  for(i in 1:length(mids)){
   p_e = dpois( mids[i],lambda=pars[1] )
   F_i = ni[i]

#   print( paste ( " F_e " ,p_e*nM , " F_i " , F_i , " (F_i - p_e*nM)^2 / (p_e*nM) " ,  (F_i - p_e*nM)^2 / (p_e*nM) ) )

   theoretical_count[i] = p_e*nM

   if( theoretical_count[i] >= 5 #nM*significance 
       ) {
     chisqr=chisqr+(F_i - p_e*nM)^2 / (p_e*nM)
     n_bins=n_bins+1
     points(mids[i],theoretical_count[i])
   } else {
     points(mids[i],theoretical_count[i],pch='x')
   }
  }


  return(list(chisqr=chisqr,n_bins=n_bins))

}



chisqr_test_pois <- function (data,lambda) {


  max_data<-max(data)

  breaks= -0.5+c(0:(max_data+2))

  hist_data <- hist(data,breaks=breaks)

#  print( paste ( " df = " , length(hist_data$counts) - 2 ))
#  print(hist_data$mids)

#  return ( chisqr_norm_2(c(mu,sigma) , hist_data$counts ,hist_data$mids ) )




  return ( chisqr_pois(c(lambda), hist_data$counts ,hist_data$mids ) )




}


chisqr_pois_sim<-function(N,k,R,lambda=2.0){

  n_breaks=round(k/lambda)

  real_lambda=k/n_breaks

  breaks=seq(0,N,length.out=n_breaks)

  chisqrs<-numeric(R)

  for( i in 1:R ) {
    sample_i<-sample(N,k)
    hist_pri<-hist(sample_i,breaks=breaks)  
    chisqrs[i] <- chisqr_test_pois(hist_pri$counts,real_lambda)
  }

  return(chisqrs)

}

