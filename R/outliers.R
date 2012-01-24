outliers <- function(data , block_len = 10 , tol_cluster = 3.0,tol_detect=5.0) {

  # length of incoming data
  len=length(data)

  ## detect where data is has nans
  nans <- is.nan(data)

  # number of block the data is grouped in
  nblocks = floor( len/block_len )

  # arrays containing the mean, standard deviation and cluster membership
  # for each block
  samples_mean <- numeric(nblocks) 
  samples_sd <- numeric(nblocks)
  samples_cluster <- rep(-1,nblocks)
  samples_skip <- rep(FALSE,nblocks)

  # store indices of global array belong to the respective block
  # this is somehow redundant but makes things easier in later stages
  indices=array(dim=c(nblocks,block_len) )
  for( i in 1:nblocks ) {
    indices[i,] = (1+block_len*(i-1)):(block_len*(i))
    if(sum(nans[indices[i,]]) > 0 ) samples_skip[i] <- TRUE
  }


  # calculate mean and sd for each block
  for( i in 1:nblocks ) {
     samples_mean[i] =  mean(data[indices[i,]])
     samples_sd[i] =    sd(data[indices[i,]])
  }
  

  # assign a cluster to each block
  cluster = 0
  for( i in 1:(nblocks-1) ) {
    
    if(samples_skip[i])  next
     
    if(samples_cluster[i] != -1 ) next
    
    cluster =  cluster + 1
    
    for( j in (i):nblocks ) {
      if( samples_skip[j] ) next
      diff = abs(samples_mean[j]-samples_mean[i] )
      if( diff < tol_cluster *samples_sd[i] && diff < tol_cluster*samples_sd[j] )
        samples_cluster[j]=cluster 
    }
  
  }

  num_clusters = cluster

  # count frequencies for each cluster
  freqs=numeric(num_clusters)
  for( i in 1:num_clusters)
    freqs[i] = sum(samples_cluster[! samples_skip]== i )
  

  # output results
##   print(paste(num_clusters , "clusters found, printing frequencies:"))
##   print(freqs)
  
  # select cluster with most members
  major_cluster <- which( freqs == max(freqs) )

  # print cluster with most members
##  print(paste(" cluster with most members is cluster " ,major_cluster ))


  # OK now we locked down
  # we know the mean and sd of the major cluster
  # this gives us a definite criterion for deciding which sample
  # "belongs" to the ensemble
  
  # calculating mean and sd of major cluster
   major_indices = which( samples_cluster == major_cluster )
   major_mean = mean( samples_mean[major_indices] )
   major_sd = mean( samples_sd[major_indices] )
##   print(paste("major_sd" , major_sd , " major_mean" , major_mean))

  # ok and now go out and find the outliers
  outlier_map = rep(FALSE,len)

  for( i in 1:len ) {
     diff = abs(data[i]-major_mean)
     if( nans [i]) {
       outlier_map[i] = TRUE
     } else  if( diff > tol_detect * major_sd) {
       outlier_map[i] = TRUE
     }
  }


  return(outlier_map)

}
