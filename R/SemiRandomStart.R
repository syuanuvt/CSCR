### obtain a semi-random cluster assignment (with a known partition)
SemiRandomStart <- function(start_part, num_exchange = 0){


  ncluster <- length(unique(start_part))
  nobs <- length(start_part)
  final_part <- start_part
  index <- 1:nobs
  if(num_exchange != 0 ){
    exchange.index <- sample(index, num_exchange, replace = FALSE)
  }
  if(num_exchange != 0 ){
    for (i in 1:num_exchange){
      if(ncluster != 2){
        ind <- exchange.index[i]
        current_cluster <- start_part[ind]
        new_cluster <- sample(setdiff(1:ncluster,current_cluster),1)
        final_part[ind] <- new_cluster
      }
      if(ncluster == 2){
        ind <- exchange.index[i]
        current_cluster <- start_part[ind]
        new_cluster <- setdiff(1:ncluster,current_cluster)
        final_part[ind] <- new_cluster
      }
    }
  }

  return (final_part)
}
