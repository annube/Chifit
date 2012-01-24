
normal_boot_error <- function(boot_result,index = 1 ) {
        boot_ci <- boot.ci(boot_result,type="norm",index=index)
        boot_error = (boot_ci$norm[3]-boot_ci$norm[2])/2./1.96
        return (boot_error);
}


normal_boot_error_from_data <- function(data,estimate){
  b = list()
b$t=array(data,dim=c(length(data),1))
b$R=length(data)
b$t0=estimate
b$call="boot"
attr(b,"class") <- "boot"
return(normal_boot_error(b))

}
