



chisqr_genral <- function (hist,f){


  n.breaks = length(hist$breaks)
  
  theoretical.counts <- numeric(length(hist$counts))

  theoretical.counts = ( sapply( hist$breaks[2:n.breaks] , function (x) f(x)) - 
    sapply( hist$breaks[1:(n.breaks-1)] , function (x) f(x)) ) * sum(hist$counts)


  range = theoretical.counts >= 5
  
  points(hist$mids[range],theoretical.counts[range])
  points(hist$mids[!range],theoretical.counts[!range],pch=4)

  chisqr = sum( ( ( hist$counts - theoretical.counts ) ^ 2 / theoretical.counts )[range] )

  return(list(chisqr=chisqr,n.intervals=sum(range)))

}
