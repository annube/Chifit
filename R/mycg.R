

mycg <- function(par,x,y,mycg.f,mycg.df,W, tol = 1.e-10,  maxsteps=50,...){
  beta=par

  rel.change=1

  steps=0

  betahist <- array(dim=c(maxsteps,length(par) ) )

  betahist[1,] = par
    
  while( rel.change>1.e-10 && steps < maxsteps){

    delta.y=y-mycg.f(x,beta,...)
    J=mycg.df(x,beta,...)

    Cw <- t(J) %*% W %*% J
    rhs <- t(J) %*% W %*% delta.y

    delta.beta <- try( solve(Cw,rhs),TRUE )

    ## do levenberg marquardt iteration if system is singular
    lambda=0.1*sum(diag(Cw))
    while( inherits( delta.beta,"try-error" ) ) {
      print(lambda)
      Cw.shifted <- Cw + diag(rep(lambda,length(beta)))
      delta.beta <- try( solve(Cw.shifted,rhs),TRUE )
      lambda=lambda*2
    }
      
    
    beta=beta+delta.beta

    rel.change = sqrt(sum(delta.beta^2)) / sqrt(sum(beta^2))
##    print(rel.change)
    steps=steps+1
    betahist[steps,]=beta
    
  }

  
    return( list( beta=beta, num.it = steps, final.prec = rel.change ,betahist=betahist[1:(steps+1),]) )

}
