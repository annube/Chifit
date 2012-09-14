simple.optim <- function(x,y,dy,f,df,par0,fit.tol = 1.e-10,maxit=100,use.line.search=FALSE,use.quadratic.regression=FALSE,upd.range=NULL,verbose=FALSE) {
  par.lm <- par0
  residue <- y-f(x,par.lm)
  Chisqr.new = sum(residue^2/dy^2)
  rel.diff = 1

  if(missing(upd.range) )
    upd.range <- 1:length(par0)
  
  it = 0
  while( abs(rel.diff)>fit.tol && it < maxit ){
    if(verbose) print(paste(Chisqr.new," ",(rel.diff)))
    dff <- df(x,par.lm)
    lm.res <- lm( residue~dff[,upd.range]+0,weights=1/dy^2)
    dpar.lm = lm.res$coefficients

    dpar.lm[ is.na(dpar.lm) ] = 0

    if( use.line.search ){
      f.line <- function(alpha) { par.test = par.lm ; par.test[upd.range] = par.lm[upd.range] + alpha*dpar.lm ; return( sum( ( ( y-f(x,par.test) )/dy ) ^2 ) ) }
      alpha.optim <- line.search(f.line,1,0.1,1.e-3,2)
      if(verbose) print(alpha.optim)
      par.lm[upd.range] = par.lm[upd.range] + dpar.lm*alpha.optim
    } else {

      if(use.quadratic.regression){
        alpha.test=c(0,0.6,1)
        
        par.new = par.lm
        par.new[upd.range] = par.lm[upd.range] + dpar.lm*alpha.test[3]
        
        par.intermed = par.lm
        par.intermed[upd.range] = par.lm[upd.range] + dpar.lm*alpha.test[2]
        
        fx=c(
          Chisqr.new,
          sum( ( ( y-f(x,par.intermed) )/dy ) ^2 ),
          sum( ( ( y-f(x,par.new) )/dy ) ^2 )
          )
        
        alpha.optim <- x0.quad.regr(alpha.test,fx)
        par.optim = par.lm
        par.optim[upd.range] = par.lm[upd.range] + dpar.lm*alpha.optim
        
        falpha.optim <- sum( ( ( y-f(x,par.optim) )/dy ) ^2 )
        
        optim.index <- order(c(fx,falpha.optim))[1]
        if(verbose) print(optim.index)
        par.lm[upd.range] = par.lm[upd.range] + dpar.lm*( c(alpha.test,alpha.optim)[optim.index] )
      } else {
        par.lm[upd.range] = par.lm[upd.range] + dpar.lm
      }
      
      
    }
    
    residue <- y-f(x,par.lm)
    Chisqr.old = Chisqr.new
    Chisqr.new = sum(residue^2/dy^2)
    rel.diff = 1 - Chisqr.new/Chisqr.old 
    it = it + 1 
  }

  if( it >= maxit ) { print("Error maxit reached " ) ; return( NA )}
  
  return(list(beta=par.lm,Chisqr=Chisqr.new,iterations=it))
}
