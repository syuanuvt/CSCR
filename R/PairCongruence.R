paircongruence <- function(lista, listb, cluster_num){

  cong <- vector("numeric", length = cluster_num)
  for (i in 1:cluster_num){
    p <- 0
    cong_i <- vector("numeric", length = (cluster_num - 1))
    current <- lista[[i]]
    for (j in 1:cluster_num){
      p <- p + 1
      cong_i[p] <- tuckercongruence(current, listb[[j]])
    }
    cong[i] <- max(cong_i)
  }
  return (cong)
}
