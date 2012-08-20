

plot.wlm <- function(wlm.res,x.index=2,ylim=NULL,...){

  l.mean <- wlm.res$predict

  l.max <- wlm.res$predict+wlm.res$dpredict
  l.min <- wlm.res$predict-wlm.res$dpredict

  if(missing(ylim))
    ylim=c( min( l.min ), max(l.max) )

  plotwitherror(wlm.res$x[,x.index],wlm.res$y,wlm.res$dy,ylim=ylim,...)
  
  lines(wlm.res$x[,x.index],l.mean)
  lines(wlm.res$x[,x.index],l.min,lty="dashed")
  lines(wlm.res$x[,x.index],l.max,lty="dashed")


  x.min = min(wlm.res$x[,x.index])
  x.max = max(wlm.res$x[,x.index])
  text(0.5 * (x.min + x.max), min((wlm.res$y-wlm.res$dy)), 
       bquote(paste(chi^2 == .(wlm.res$Chisqr), "; ", dof == .(wlm.res$dof))))
  
}
