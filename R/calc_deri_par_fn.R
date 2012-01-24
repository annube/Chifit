

convert.expr.vars.fw <- function( expr.char , var.name ){

  uid = "VpG4nWJvlQ57t4Nt190O"

  var.name = gsub( "[.]","[.]" , var.name )

  brace.open.pat = "[[:space:]]*\\[[[:space:]]*"
  brace.close.pat = "[[:space:]]*\\]"
  comma.pat = "[[:space:]]*,[[:space:]]*"

  ## one index
  parpattern1 = sprintf( "(%s)%s([[:alnum:]]*)%s" , var.name , brace.open.pat,brace.close.pat )
  expr.char.mod = gsub(parpattern1,paste("\\1.",uid,".\\2",sep=""),expr.char)

  ## two indices
  parpattern2 = sprintf( "(%s)%s([[:alnum:]]*)%s([[:alnum:]]*)%s" , var.name ,brace.open.pat,comma.pat,brace.close.pat )
  expr.char.mod = gsub(parpattern2,paste("\\1.",uid,".\\2_c_\\3",sep=""),expr.char.mod)

  ## three indices
  parpattern3 = sprintf( "(%s)%s([[:alnum:]]*)%s([[:alnum:]]*)%s([[:alnum:]]*)%s" , var.name ,brace.open.pat,comma.pat,comma.pat,brace.close.pat )
  expr.char.mod = gsub(parpattern3,paste("\\1.",uid,".\\2_c_\\3_c_\\4",sep=""),expr.char.mod)
  
  return(expr.char.mod)
}

convert.expr.vars.bw <- function( expr.char.mod , var.name ){
  uid = "VpG4nWJvlQ57t4Nt190O"

  var.name = gsub( "[.]","[.]" , var.name )

  s.pattern3 = sprintf( "(%s).%s.([[:alnum:]]*)_c_([[:alnum:]]*)_c_([[:alnum:]]*)" , var.name , uid )
  expr.char = gsub(s.pattern3,paste("\\1[\\2,\\3,\\4]",sep=""),expr.char.mod) 

  s.pattern2 = sprintf( "(%s).%s.([[:alnum:]]*)_c_([[:alnum:]]*)" , var.name , uid )
  expr.char = gsub(s.pattern2,paste("\\1[\\2,\\3]",sep=""),expr.char) 

  s.pattern1 = sprintf( "(%s).%s.([[:alnum:]]*)" , var.name , uid )
  expr.char = gsub(s.pattern1,paste("\\1[\\2]",sep=""),expr.char) 
  
  
  return( expr.char )
  
}


## calculates the derivative of a parametric function
## given by expr with respect to the parameters
calc.deri.par.fn <- function( expr , parlen, additional.pars = c("x") ) {
  uid = "VpG4nWJvlQ57t4Nt190O"

  names = sapply(1:parlen,function(x) sprintf("par.%s.%d",uid,x) )

  names.stub = sprintf("par.%s.",uid)

  
  expr.char = deparse( expr,width.cutoff = 500 )

  parpattern = "(par)[[:space:]]*\\[[[:space:]]*([0-9]+)[[:space:]]*\\]"
  
  expr.char.mod = gsub(parpattern,paste("\\1.",uid,".\\2",sep=""),expr.char)


  expr.mod = parse(text = expr.char.mod)
  


  deris <- sapply(names,function(x) D(expr.mod,x) )

  deris.deparse <- sapply(deris,function(x) deparse(x,width.cutoff=500))

  
  deris.deparse.mod <-
    sapply( deris.deparse , function(x) 
           gsub( paste(names.stub , "([0-9]+)" ,sep = "" ) , "par\\[\\1\\]", x )
           )

  
  names(deris.deparse.mod) =  sapply(1:parlen,function(x) sprintf("par[%d]",x) )
  
  return( parse( text = deris.deparse.mod )  )
}


test.calc.deri.par.fn <- function(f,parlen) {
  ##  f <-  function(par,x) par[1] * exp( - par[2] * x )
  
  expr <- if (identical(body(f)[[1]], as.name("{"))) body(f)[[2]] else body(f)
  
  df.exp <- calc.deri.par.fn( expr , parlen )
  df <- function(par,x) sapply( df.exp , function(y ) eval(y,list(par=par,x=x)) )
  
  return(df)
}


