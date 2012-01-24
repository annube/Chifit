



make.level.indices <- function(array){
  levels <- as.numeric(levels(as.factor(array)))
  indices <- numeric(length(array))

  for( i in 1:length(levels))
    indices[ which(array == levels[i]) ] = i
  
  return(indices)
}
  
