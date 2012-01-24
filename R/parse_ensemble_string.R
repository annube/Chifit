parse_ensemble_string <- function(e_string){
  
  ssplit <- unlist(strsplit(e_string,'[.]'))
  ssplit1 <- unlist(strsplit(ssplit[1],NULL))
  beta <- c(1.9,1.95,2.1) [
                           which( c("A","B","D") == ssplit1[1] ) 
                           ]
  mu <- as.numeric(paste(ssplit1[2:length(ssplit1)],collapse=""))/10000

  
  search.res <- regexpr("[0-9]{2,3}",ssplit[2])
  matchlen <- attr(search.res,"match.length")
##  print(matchlen)

  ssplit2 <-  unlist(
                     strsplit(ssplit[2],NULL)
                   )

  
  L= as.integer(
    paste(
          ssplit2
                   [search.res:(matchlen)] ,collapse=""
          )
    )
  

  if( (matchlen+1) <= length(ssplit2) ) {
    modifier= 
      paste(
            ssplit2
            [(matchlen+1):length(ssplit2)] ,collapse=""
            )
  } else {
    modifier <- "-"
  }

  
##  print( noquote(paste("beta = ",beta," mu = " ,mu , " L= " , L ," mod = \"" , modifier,'"',sep="")))

  return(list(beta=beta,mu=mu, L=L, modifier = modifier))
}
