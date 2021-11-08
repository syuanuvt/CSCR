### obtain a random cluster assignment (typically as the first step)
RandomStart <- function(mem_cluster){

  i <- sum(mem_cluster)
  ncluster <- length(mem_cluster)

  mem.cluster <- matrix(0, nrow = i, ncol = ncluster)
  ind.cluster <- vector("numeric")

  # initialize
  for (t in 1:ncluster){
    ind.cluster <- c(ind.cluster, rep(t, mem_cluster[t]))
  }

  # randomlize
  ind.cluster <- sample(ind.cluster)

  # set to the membership cluster
  for (j in 1:i){
    mem.cluster[j, ind.cluster[j]] <- 1
  }

  return (list(mem = mem.cluster, ind = ind.cluster))
}
