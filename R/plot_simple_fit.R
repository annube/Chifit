

## library(hadron)

plot_simple_fit <- function(data_x,data_y,data_dy,fitfn,file="simple_fit.pdf")
{

  pdf(file=file)
  plotwitherror(data_x,data_y,data_dy)

  xs = seq(min(data_x),max(data_x),length.out=100)
  ys = fitfn(xs)

  lines(xs,ys)
  dev.off()
  
}

