#' Create all possible partitions
#'
#' This function creates all possible partitions of component structure
#'
#' @param com_total The total number of components
#' @param nblock The total number of data blocks
#' @return A matrix where each column indicates a specific combination of components
#' @export
#' 
AllPartition <- function(com_total, nblock){
  partition <- list()
  num <- 0
  comb <- as.data.frame(permutations((com_total+1),nblock,0:com_total, set = FALSE, repeats.allowed = TRUE))
  comb <- cbind(comb, (com_total - apply(comb[,1:nblock],1,sum)))
  comb_left <- comb[comb[,(nblock+1)]>=0,]
  names(comb_left)[ncol(comb_left)] <- "last"
  return(comb_left)
}