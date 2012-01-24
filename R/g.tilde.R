
g.tilde.c <- function(x,n=1){
	.Call( "g_tilde",x,n, PACKAGE = "chifit" )
}
