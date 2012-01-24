




simulate_fit <- function( parameter,fn,input_vals ,N=200 , ... )
{



n_params=length(parameter)


pars = array( N * n_params, dim = c (n_params,N) )

chisqrs = numeric(N)


variances = input_vals$data_dy # /1.96
means = input_vals$data_y

  sim_input=input_vals



for ( i in 1:N ) {

  for( j in 1:n_params ) sim_input$data_y[j] = rnorm(1,mean=means[j],sd=variances[j])

  sim_optim_fit <- optim(par=parameter,
                	fn=fn,
                 	input_vals=sim_input,
                        ...
                	)

  pars[,i] = sim_optim_fit$par

  chisqrs[i] = sim_optim_fit$value

  if( i %% (N/4) == 0 ) print(sprintf(" %5.1f%% finished", i/N*100))

}

  optim_fit <- optim(par=parameter,
                	fn=fn,
                 	input_vals=input_vals,
                        ...
                	)



  return(list(optim_fit=optim_fit,simulation=pars,chisqrs=chisqrs,fit_fn=fn,input_vals=input_vals))


}


