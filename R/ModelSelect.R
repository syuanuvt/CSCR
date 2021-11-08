#### select the optimal number of components and clusters
ModelSelect <- function(data_x, data_y, nvar,nblock, r_max, k_max){
  ########### select the total number of components
  results <- list()
  res_x <- rep(0,r_max)
  for (i in 1 : r_max){
    results[[i]] <- csca_cpp(data_x,nvar, nblock, i, 1, 200)
    res_x[i] <- results[[i]]$loss
  }
  opt_results <- ScreeSelect(res_x)
  opt_com <- opt_results[[1]]
  ########## select the total number of clusters
  results_loss <- list()
  loss_x <- rep(0,k_max)
  for (k in 1:k_max){
    results <- csca_cpp(data_x,nvar, nblock, opt_com, k, 200)
    results_loss[[k]] <- MultiCSCCR(data_x,data_y,nblock, opt_com, c(0,0), nvar, k, alpha = .99,
                                    converge = 1e-3, iteration = 100, num_starts = 5, type = "known", con = TRUE,
                                    start_part = results$cluster_mem)
    loss_x[k] <- results_loss[[k]]$loss
  }
  opt_results_cluster <- ScreeSelect(loss_x)
  opt_cluster <- opt_results_cluster[[1]]
  ######## go back to select the number of components
  results <- list()
  res_x <- rep(0,r_max)
  for (i in 1 : r_max){
    results[[i]] <- csca_cpp(data_x,nvar, nblock, i, opt_cluster, 200)
    res_x[i] <- results[[i]]$loss
  }
  opt_results <- ScreeSelect(res_x)
  opt_com <- opt_results[[1]]
  #########
  all_comb <- AllPartition(opt_com, nblock)
  loss_str <- rep(0, nrow(all_comb))
  index_str <- rep(0, nrow(all_comb))
  for (i in 1:nrow(all_comb)){
    n_com <- all_comb[i,1]
    n_distinct <- all_comb[i,-1]
    if(i == nrow(all_comb)){
      results_str[[i]] <- csca_cpp(data_x,nvar, nblock, n_com, opt_cluster, 200)
    }
    if(i != length(ncom_all)){
      results_str[[i]] <- cssca_quick_cpp(data_x, nvar, nblock, n_com, n_distinct, opt_cluster,
                                          nrespondents = nrow(data_x), sparsity = 0, feed = results_loss[[opt_cluster]]$cluster_assign, 1/6)
    }
    loss_str[i] <- results_str[[i]]$loss
    index_str[i] <- all_var - loss_str[i]
  }

  selection <- list(opt_com = opt_com, opt_cluster = opt_cluster, loss_str = loss_str, index_str = index_str)
  return(selection)
}
