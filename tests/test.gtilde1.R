library(chifit)
.Call("eval_gtilde1_R",2.,PACKAGE="chifit")
g.tilde(2.,partitions)

.Call("eval_dgtilde1_R",2.,PACKAGE="chifit")
g.tilde.deri(2.,partitions)
