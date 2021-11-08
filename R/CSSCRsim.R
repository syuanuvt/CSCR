CSSCRSim <- function(ncluster, memcluster, ncom, ndistinct, nblock, nvar, ynvar,
                     psparse = 0, pnoise_x = 0, pnoise_y = 0, clust_reg, com_rele, comp_str = NULL,
                     equal_loading = TRUE, n_test = NULL,
                     pmean_x = 0, pmean_y = 0, p_combase = 0) {

  # ncluster: the number of cluster
  # memcluster: number of participants in each cluster
  # ncom: number of common components
  # ndistinct: number of distinct components (a vector of length nblock, with the ith element indicating the number of distinctive components
  # in ith data block)
  # nblcok: number of blocks
  # nvar: number of variables (vector indicates the number of variables for in each data block)
  # ynvar: the number of variables in the outcome block
  # psparse: percentage of sparsity (vector indicates the percentage for each block)
  # pnoise_x: percentage of noise for predictors
  # pnoise_y: percentage of noise for the outcome
  # clust_reg: cluster-specific regression coefficients
  # com_rele: the relevance of the components (the length equals the total number of components)
  # n_test: number of tests per data block


  # irregular input
  if (length(memcluster) == 1){
    memcluster <- rep(memcluster, ncluster)
  }

  if (length(nvar) == 1){
    nvar <- rep(nvar, nblock)
  }

  if (length(ndistinct) == 1 & sum(ndistinct) != 0){
    ndistinct <- rep(ndistinct, nblock)
  }

  # the vector indicates the relationship between distinctive component and the number of data block
  distinct_index <- vector("numeric", length = sum(ndistinct))
  sum_var <- 0
  ini_var <- 0
  common <- ncom
  pos <- 1
  for (i in 1:nblock){
    if (sum(ndistinct) != 0){
      distinct_index[pos: (pos + ndistinct[i] - 1)] <- i
      pos <- pos + ndistinct[i]
    }
    ini_var <- ini_var + nvar[i]
    sum_var <- c(sum_var, ini_var)
  }

  # the aggregate level of information
  all_component <- sum(common, ndistinct)
  all_member <- sum(memcluster)
  all_var <- sum(nvar)
  if(is.null(comp_str)){
    comp_str <- rep(1/all_component,all_component)
  }

  # for the easy computation (fill in the 0th item)
  cluster_rep <- c(0, memcluster)
  block_rep <- c(0, nvar)

  # cluster assignment (becuase of the randomness, the assignment could be arbitrary, just from the first avilable
  # observation looping towards the last avilable observation)
  cluster_mem <- vector("numeric", length = all_member)
  for  (i in 1:ncluster){
    cluster_mem[(sum(cluster_rep[1:i]) + 1):sum(cluster_rep[1:i+1])] <- i
  }

  ## create the initial score matrix
  true_score <- list()
  for (k in 1:ncluster){
    indices <- which(cluster_mem == k)
    ori_score <- matrix(rnorm(length(indices) * all_component), nrow = length(indices), ncol = all_component)
    score_processed <- MatrixCenter(ori_score, 1, 1)
    score.de <- svd(score_processed)
    true_score[[k]] <-  score.de$u %*% t(score.de$v) * sqrt(nrow(ori_score)) #the alternative way to create the starting score matrix (that are orthogonal)
  }
  # the loading matrix for the predictors

  px <- list()
  if (equal_loading == FALSE){
    com_base <- sqrt(p_combase) * matrix(runif(all_var * all_component, min = -1, max = 1), nrow = all_var, ncol = all_component)
    for (i in 1:ncluster){
      px[[i]] <- com_base + sqrt(1 - p_combase) * matrix(runif(all_var * all_component, min = -1, max = 1), nrow = all_var, ncol = all_component)
    }
  }
  if(equal_loading == TRUE){
    px_all <- matrix(runif(all_var * all_component, min = -1, max = 1), nrow = all_var, ncol = all_component)
    for (i in 1:ncluster){
      px[[i]] <- px_all
    }
  }

  # generate the positions of structure-induced zeros
  distinct_zeros <- vector("numeric")
  if (sum(ndistinct) != 0){
    num_var <- 0
    for (i in 1:sum(ndistinct)){
      distinct_zeros <- c(distinct_zeros, ((all_var * (common + i - 1)) + sum_var[distinct_index[i]] + 1): ((all_var * (common + i - 1)) + sum_var[(distinct_index[i] + 1)]))
    }
  }
  for (i in 1:ncluster){
    px[[i]][distinct_zeros] <- 0
    px[[i]] <- px[[i]] * sqrt(1-pmean_x) * sqrt(1-pnoise_x)/sqrt(mean(apply(px[[i]]^2,1,sum)))
    each_variance <- apply(px[[i]]^2,2,sum)
    all_variance <- sum(each_variance)
    for(j in 1:all_component){
      px[[i]][,j] <- px[[i]][,j]/sqrt(each_variance[j])*sqrt(all_variance*comp_str[j])
    }
  }

  # true score for the regression problem
  #true_score_y <- true_score)
  final_x <- matrix(nrow = all_member, ncol = all_var)
  final_y <- matrix(nrow = all_member, ncol = ynvar)
  py <- list()

  for (k in 1:ncluster){
    py[[k]] <- matrix(nrow = ynvar, ncol = all_component)
    for(i in 1:all_component){
      py[[k]][1,i] <- clust_reg[k,i]*com_rele[i]
    }
    if(ynvar > 1){
      for (i in 2:ynvar){
        for(j in 1:all_component)
          py[[k]][i, j] <- rnorm(1,py[[k]][1, j],1)
      }
    }
  }

  true_score_all <- matrix(nrow = all_member, ncol =all_component)
  for (k in 1:ncluster){
    indices <- which(cluster_mem == k)
    true_score_k <- true_score[[k]]
    true_score_all[indices, ] <- true_score_k

    ## create the final data
    final_k <- true_score_k %*% t(px[[k]])
    final_k <- MatrixCenter(final_k,1,0)
    final_x[indices, ] <- final_k

    final_k <- true_score_k %*% t(py[[k]])
    final_y[indices,] <- final_k- mean(final_k)
  }

  test.index <- c()
  if(!is.null(n_test)){
    for (i in 1:ncluster){
      test.index <- c(test.index, ((i-1)*memcluster[1]+1):((i-1)*memcluster[1]+n_test))
    }
  }
  train.index <- setdiff(1:(all_member),test.index)
  #final_y <- final_y * sqrt(1-pmean_y) * sqrt(1-pnoise_y)/mean(final_y^2)

  add_noise_x <- Add(final_x, cluster_mem, pnoise_x, pmean_x)
  add_noise_y <- Add(matrix(final_y, ncol = 1), cluster_mem, pnoise_y, pmean_y)
  noise_x <- add_noise_x[[2]]
  noise_y <- add_noise_y[[2]]
  final_x <- add_noise_x[[1]]
  final_y <- add_noise_y[[1]]
  mean_y <- add_noise_y[[3]]
  mean_x <- add_noise_x[[3]]
  final_x_test <- final_x[test.index,]
  final_x <- final_x[-test.index,]
  final_y_test <- final_y[test.index]
  final_y <- final_y[-test.index]

  sim = list(x = final_x, y = final_y, cluster_mem = cluster_mem[train.index], score_all = true_score_all[train.index,],
             loading_x = px, loading_y = py, noise_x = noise_x, noise_y = noise_y,
             x_test = final_x_test, y_test = final_y_test, cluster_mem_test = cluster_mem[test.index],
             mean_y = mean_y, mean_x = mean_x)
  return(sim)
}
