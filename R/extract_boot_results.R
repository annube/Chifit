extract_boot_results <- function(index_table,index=1){
num_masses = length(index_table$i)
data_x = numeric(num_masses)
data_x0 = numeric(num_masses)
data_y = numeric(num_masses)
data_dy = numeric(num_masses)
num_boot_samples = NA
data_y_boot = NA
first = 1

## construct a table of results from the boot samples
for( i in 1:num_masses ) {
  
  boot_quantities_file<-sprintf("boot_quantities.%02d.%02d.Rdata",index_table$i[i],index_table$j[i])
  load(boot_quantities_file)

 if( first == 1) {
    first = 0
    num_boot_samples = length(boot_quantities$t[,index])
    data_y_boot = array(NA,dim=c(num_boot_samples,num_masses))
  }

  data_x[i] = masses$V1[index_table$j[i]+1]
  data_x0[i] = masses$V1[index_table$i[i]+1]
  data_y[i] = boot_quantities$t0[index] # index = 1 -->>  m_eff
  data_dy[i] = normal_boot_error(boot_quantities,index=index)
  data_y_boot[,i] = boot_quantities$t[,index]

}
 
 return(list(data_x=data_x,data_x0=data_x0,data_y=data_y,data_dy=data_dy,data_y_boot=data_y_boot))
 
}