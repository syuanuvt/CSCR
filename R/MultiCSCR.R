MultiCSCR <- function(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster, alpha = .8,
                       converge = 1e-8, iteration = 1000, num_starts = 1, type = NULL, con = FALSE,
                       start_part = NULL){

  nobs <- nrow(data_x)
  start_partition <- list()
  if(n_cluster !=1 ){
    mem_cluster <- rep(nrow(data_x)/n_cluster, n_cluster)
    if(type == "random"){
      for(i in 1:num_starts){
        rs.results <- RandomStart(mem_cluster)
        start_partition[[i]] <- rs.results$ind
      }
    }
    if(type == "known"){
      num <- floor(nobs/10)
      for(i in 1:num_starts){
        if(i == 1){
          start_partition[[i]] <- start_part
        }
        if(i != 1){
          start_partition[[i]] <- SemiRandomStart(start_part, num)
        }
      }
    }
  }
  if(n_cluster == 1){
    start_partition[[1]] <- start_part
    num_starts <- 1
  }
  for (i in 1:num_starts){
    if(i == 1){
      if(con == FALSE){
        est <- CSSCR_nocon(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster,alpha,
                           converge, converge, iteration, start_part = start_partition[[i]])
        final_loss <- est$loss
        final_results <- est
      }
      if(con == TRUE){
        est <- CSSCR_group(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster, alpha,
                           converge, converge,  iteration, start_part = start_partition[[i]])
        final_loss <- est$loss
        final_results <- est
      }
    }
    if(i != 1){
      if(con == FALSE){
        est <- CSSCR_nocon(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster, alpha,
                           converge, converge, iteration, start_part = start_partition[[i]])
        current_loss <- est$loss
        if(current_loss < final_loss){
          final_results <- est
          final_loss <- current_loss
        }
      }
      if(con == TRUE){
        est <- CSSCR_group(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster, alpha,
                           converge, converge, iteration, start_part = start_partition[[i]])
        current_loss <- est$loss
        if(current_loss < final_loss){
          final_results <- est
          final_loss <- current_loss
        }
      }
    }
  }
  final_results$loss <- final_loss
  return(final_results)
}
