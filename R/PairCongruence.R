#' Calculate the congruence of two list of matrices
#'
#' The function calculates the congruence of two lists as indicated by the average Tucker's congruence
#'
#' @param lista the first list of all matrices
#' @param listb the second list of all metrics
#' @param cluster_num the number of clusters
#' @return A scale indicating the average congruence level of the two lissts
#' @export
#' 
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
