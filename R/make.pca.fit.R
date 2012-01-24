make.pca.fit <- function(Cc,y,dy,x.tuple,n.pca.y,n.pca.x) {

  Cc.eig <- eigen(Cc)

##  print(paste("condition number of Cc: " , max(Cc.eig$values)/min(Cc.eig$values) ))


   y.pca.all = as.vector ( y %*% Cc.eig$vectors )
   dy.pca = Cc.eig$vectors %*% diag(y.pca.all)
   r <- apply(dy.pca,2,function(x) sum((x-y)^2)/sum(y^2) )
   pca.indices <- order(r)[1:n.pca.y]


  P.y <- Cc.eig$vectors[,1:n.pca.y] #pca.indices]
#  P.y <- Cc.eig$vectors[,pca.indices]


  y.pca = y %*% P.y


  y.check = P.y %*% as.vector(y.pca)

print(  max( abs(as.vector(y.check)-y)/dy ) )

                                        # calculate reduced correlation matrix
                                        #C <- cov((t(y.boot) %*% P.y))
  C = t(P.y) %*% Cc %*% P.y

  dy.pca = sqrt(diag(C))
  C.eig <- eigen(C)
##  print(paste("condition number of C: " , max(C.eig$values)/min(C.eig$values) ))


  ## reduce to pca in x

  W <- solve(C)

  X <- t(P.y) %*% t(dmultipoly(1,x.tuple)[c(1:24,26:28),])

  Cw <- t(X) %*% W %*% X
  Cw.eig <- eigen(Cw)
##  print(paste("condition number of Cw: " , max(Cw.eig$values)/min(Cw.eig$values) ))


  P.x = Cw.eig$vectors[,1:n.pca.x]

  X.pca <- X %*% P.x

  Cw <- t(X.pca) %*% W %*% X.pca
  Cw.eig <- eigen(Cw)

##  print(paste("condition number of Cw: " , max(Cw.eig$values)/min(Cw.eig$values) ))


  ## now do the fits in the reduced space

  lm.fit <- lm(as.vector(y.pca)~X.pca+0)

  y.pca.predict <- X %*% P.x %*% lm.fit$coefficients

  y.predict <- as.vector ( P.y %*% y.pca.predict )

                                        # plotwitherror(1:length(y),as.vector(y.predict) - y,dy)

  sum( ((y.predict-y)/dy)^2 )



  wlm.res <- wlm(X.pca,as.vector(y.pca),C=C)#,dy=sqrt(diag(C)))

  wlm.res$Chisqr

  y.pca.predict <- X %*% P.x %*% wlm.res$beta

  y.predict <- as.vector(P.y %*% y.pca.predict)
  y.res=y-y.predict
  W <- diag(1/dy^2)
  lm.chisqr = sum( y.res*as.vector( W %*% y.res) )


                                        # plotwitherror(1:length(y),y.predict - y,dy)


  
  par.wlm=P.x %*% wlm.res$beta
  par.lm=P.x %*% lm.fit$coefficients


  
  return(list( wlm.res=wlm.res,
              lm.fit=lm.fit,
              par.lm=par.lm,
              par.wlm=par.wlm,
              lm.chisqr = lm.chisqr,
              P.x=P.x,P.y=P.y
              )
         )
}
