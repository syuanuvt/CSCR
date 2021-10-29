## the current script is a test of the novel CSSCR method
library(pracma)
library(remotes)
library(ggplot2)
library(GGally)
library(dplyr)
library(tidyr)
library(rARPACK)
library(RSpectra)
library(pracma)
library(psych)
library(combinat)
library(gtools)
library(mbclusterwise)
library(flexmix)
library(combinat)
library(mclust)
library(fpc)
library(iCluster)
library(ClusterSSCA)
library(PCovR)
########## several supplementary functions
### should be run before running the main function
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
###############################################
AllPartition <- function(com_total, nblock){
  partition <- list()
  num <- 0
  comb <- as.data.frame(permutations((com_total+1),nblock,0:com_total, set = FALSE, repeats.allowed = TRUE))
  comb <- cbind(comb, (com_total - apply(comb[,1:nblock],1,sum)))
  comb_left <- comb[comb[,(nblock+1)]>=0,]
  names(comb_left)[ncol(comb_left)] <- "last"
  return(comb_left)
}
###############################################
######### function to select the optimal value of the parameter with a scree ratio test
ScreeSelect <- function(loss){
  loss_len <- length(loss)
  ratio_x <- rep(0,(loss_len-2))
  for (i in 2:(loss_len-1)){
    ratio_x[(i-1)] <- (loss[i+1] - loss[i])/(loss[i] - loss[i-1])
  }
  opt <- which.min(ratio_x)+1
  opt_1 <- which(diff(diff(ratio_x))<0)[1]+1
  results <- list(opt, opt_1, ratio_x)
  return (results)
}
### centering and scaling the matrix
MatrixCenter <- function(matrix, center, scale){
  # matrix: the matrix that need to be centered and (or) scaled
  # center: does the matrix need to be centered (1 = Yes, 0 = No)
  # scaled: does the matrix need to be scaled (1 = Yes, 0 = No)

  variable_mean <- apply(matrix, 2, mean)
  variable_sd <- apply(matrix, 2, sd)
  n_observation <- nrow(matrix)

  if (center == 1){
    matrix <- matrix - rep(1, n_observation) %*% t(variable_mean)
  }
  if (scale == 1){
    matrix <- t(t(matrix) / variable_sd)
  }

  return(matrix)
}

#### congruence
tuckercongruence <- function(matrix1, matrix2){
  ## note that two matrices should be in the same size
  # matrix1: reference matrix
  # matrix2: matrix to compare to the reference

  m <- nrow(matrix1)
  n <- ncol(matrix2)

  indic_perms <- permn(1:n)
  tuck <- vector("numeric")
  for (i in 1:length(indic_perms)){
    matrix2_perm <- matrix2[ ,indic_perms[[i]]]
    tuck[i] <- tr(abs(factor.congruence(matrix1, matrix2_perm))) / n
  }

  tuck_max <- max(tuck)
  return (tuck_max)
}

### a list of congruence
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

### add the noise to the dataset
Add <- function(inidata, cluster_assign, p_noise, p_mean){
  # inidata: the data generated via the covariance structure
  # p_noise: pencentage of noise

  n_row <- nrow(inidata)
  n_col <- ncol(inidata)
  n_cluster <- length(unique(cluster_assign))

  inidata <- MatrixCenter(inidata, center =1, scale=0)
  finaldata <- matrix(nrow = n_row, ncol = n_col)
  finaldata <- inidata
  if(n_col == 1){
    var.unit <- mean(inidata^2) / ((1-p_mean) * (1 - p_noise))
  }
  
  mean_data <- matrix(nrow = n_row, ncol = n_col)
  mean_cluster <- matrix(nrow = n_cluster, ncol = n_col)
  
  for (i in 1:n_cluster){
    mean_cluster[i, ] <- rnorm(n_col)
    mean_data[which(cluster_assign == i), ] <- rep(mean_cluster[i, ], each = sum(cluster_assign == i))
  }
  mean_data <- MatrixCenter(mean_data,1,0)
  
  if(n_col != 1){
    finaldata <- inidata + t(sqrt(p_mean * (1 - p_noise) / apply(mean_data, 2, function(x) mean(x^2))) * t(mean_data))
  }
  if(n_col == 1){
    finaldata <- inidata + sqrt(p_mean * (1 - p_noise) * var.unit / mean(mean_data^2)) * mean_data
  }
  if(n_col == 1){
    final_mean_data <- sqrt(p_mean * (1 - p_noise) * var.unit / mean(mean_data^2)) * mean_data
  }
  if(n_col != 1){
    final_mean_data <- sqrt(p_mean * (1 - p_noise) / mean(mean_data^2)) * mean_data
  }

  noise <- matrix(rnorm(n_row * n_col), n_row, n_col)
  if(n_col == 1){
    noise <- sqrt(p_noise * var.unit /mean(noise^2)) * noise
  }
  if(n_col != 1){
    noise <- t(sqrt(p_noise/apply(noise,2,function(x) mean(x^2))) * t(noise))
  }
  finaldata <- finaldata + noise

  return (list(final_data = finaldata, noise = noise, final_mean_data = final_mean_data))
}

### obtain a random cluster assignment (typically as the first step)
RandomStart <- function(mem_cluster){
  #i: number of observations
  #ncluster: number of clusters
  #ncom: number of components
  #i <- 20
  #ncluster <- 4
  #set.seed(seed)

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

######calculate the real partition of the observations
RealPar <- function(n_cluster, data_y,true_score,loading_y){
  
  n_obs <- length(data_y)
  res_y <- matrix(nrow = n_obs, ncol = n_cluster)
  real_par <- rep(NA, length = n_obs)
  danger_par <- rep(0, length = n_obs)
  for (k in 1:n_cluster){
    res_y[,k] <- data_y - true_score %*% t(loading_y[[k]])
  }
  for (i in 1:n_obs){
    real_par[i] <- which.min(abs(res_y[i,]))
    #if(sum(res_y[i,]^2 > (sum(noise_y^2) / n_obs))==0){
    #  danger_par[i] <- 1
    #}
  }
  
  return(list(real_cluster_mem = real_par, ucluster_distance = res_y) )
}

ncluster <- 2
memcluster <- c(55,55)
ncom <- 1
ndistinct <- c(0,0)
nblock <- 2
nvar <- c(1,1)
ynvar <- 1
psparse <- 0
pnoise_x <- .1
pnoise_y <- 0
clust_reg <- as.matrix(c(1,-1))
com_rele <- 1
comp_str <- 1
equal_loading <- FALSE
n_test <- 10
pmean_x <- .1
pmean_y <- .1
p_combase <- 0
### data simulation with CSSCR (also include the relevance of the components)
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
  # the positions that are not yet specified as zeros
  #retain_zeros <- setdiff(1: (all_var * all_component), distinct_zeros)

  # generate the component loading matrix for the predictors
  # the number of zeros in component loadings
  #num_zeros <- round(length(retain_zeros) * psparse)

  # determine the places for sparsity-imposed zeros
  #sparse_zeros <- sample(retain_zeros, num_zeros)
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

### without any model selection
CSSCR <- function(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster, alpha = .8,
                  converge = 1e-7, converge_2 = 1e-7, iteration = 1000, start_part = NULL){
  
  ## fix the sparsity level to zero (p_sparse = 0)
  upper <- 1e9
  stop <- 0
  iter <- 0
  
  all_member <- nrow(data_x)
  sum_var <- ncol(data_x)
  y_var <- ncol(data_y)
  block_version_data <- list(length = 2)
  block_version_data[[1]] <- data_x[,1:n_var[1]]
  block_version_data[[2]] <- data_x[, (n_var[1] + 1):n_var[2]]
  n_total <- sum(n_com, n_distinct)
  
  ## structure-induced zeros
  distinct_index <- vector("numeric")
  distinct_zeros <- vector("numeric")
  if (sum(ndistinct) != 0){
    
    all_var <- 0
    ini_var <- 0
    for (p in 1:n_block){
      distinct_index <- c(distinct_index, rep(p, n_distinct[p]))
      ini_var <- ini_var + n_var[p]
      all_var <- c(all_var, ini_var)
    }
    
    for (r.distinct in 1:sum(n_distinct)){
      distinct_zeros <- c(distinct_zeros, ((sum_var * (n_com + r.distinct - 1)) + all_var[distinct_index[r.distinct]] + 1): ((sum_var * (n_com + r.distinct - 1)) + all_var[(distinct_index[r.distinct] + 1)]))
    }
  }
  
  # set the upper bound of minimum loss
  loss_all <- vector("numeric", length = iteration)
  
  # the starting partition
  start <- vector("numeric", length = n_cluster)
  for (y in 1:n_cluster){
    start[y] <- round(all_member / n_cluster)
  }
  start[n_cluster] <- all_member - (n_cluster - 1) * round(all_member / n_cluster)
  
  
  py <- list()
  if(is.null(start_part)){
    cluster_assign <- RandomStart(start)[[2]]
  }
  if(!is.null(start_part)){
    cluster_assign <- start_part
  }
  
  ###################################################
  #cluster_assign <- a
  svd_data <- svds(data_x, n_total)
  t <- svd_data$u
  x_square <- sum(data_x^2)
  y_square <- sum(data_y^2)
  #beta <- .5
  #alpha <- .5
  beta <- alpha * y_square / (alpha * y_square + (1-alpha) * x_square)
  
  #data_x <- MatrixCenter(data_x,1,0)
  #beta <- 0.9
  
  loss_min <- upper
  loss_p_all <- rep(NA, iteration)
  loss_t_all <- rep(NA, iteration)
  loss_t_all_1 <- rep(NA, iteration)
  loss_t_all_2 <- rep(NA, iteration)
  loss_t_all_3 <- rep(NA, iteration)
  for (v in 1:iteration){
    if(loss_min < converge & v > 20)  break
    loss_new <- rep(NA, all_member)
    if (v == 1){
      px <- t(data_x) %*% t
      px[distinct_zeros] <- 0
    }
    beta * sum((data_x - t %*% t(px))^2)
    #index.nsmallest <- sort(px^2, index.return = TRUE)$ix[1:(sparse_zeros+length(distinct_zeros))]
    #x[index.nsmallest] <- 0
    
    if (v == 1){
      for(k in 1:n_cluster){
        cluster_k <- which(cluster_assign == k)
        t_k <- t[cluster_k,]
        y_k <- data_y[cluster_k,]
        py[[k]] <- t(y_k) %*% t_k %*% inv(t(t_k) %*% t_k)
      }
    }
    
    ############################
    
    t <- matrix(0, nrow = all_member, ncol = n_total)
    for (k in 1:n_cluster){
      cluster_k <- which(cluster_assign == k)
      x_k <- data_x[cluster_k,]
      xp <- beta * (x_k %*% px)
      y_k <- data_y[cluster_k,]# - mean(data_y[cluster_k])
      xp1 <- (1-beta) * (y_k %*% py[[k]])
      xp_update <- xp + xp1
      xp_svd <- svd(xp_update)
      t[cluster_k,] <- (xp_svd$u %*% t(xp_svd$v)) / sqrt(n_cluster)
    }
    ####################################
    loss_t2 <- 0
    for(k in 1:n_cluster){
      cluster_k <- which(cluster_assign == k)
      t_k <- t[cluster_k,]
      y_k <- data_y[cluster_k,]#-mean(data_y[cluster_k])
      loss_k <- (1-beta) * sum((y_k - t_k %*% t(py[[k]]))^2)
      loss_t2 <- loss_t2+loss_k
    }
    
    px <- t(data_x) %*% t
    px[distinct_zeros] <- 0
    
    ######################################
    loss_t <- beta * sum((data_x - t %*% t(px))^2)
    #loss_t_all_1[v] <- loss_t
    loss_t_all_2[v] <- 0
    for(k in 1:n_cluster){
      cluster_k <- which(cluster_assign == k)
      t_k <- t[cluster_k,]
      y_k <- data_y[cluster_k,]#-mean(data_y[cluster_k])
      reg_results <- t(y_k) %*% t_k %*% inv(t(t_k) %*% t_k)
      py[[k]] <- reg_results
      loss_k <- (1-beta) * sum((y_k - t_k %*% t(py[[k]]))^2)
      loss_t <- loss_t+loss_k
      loss_t_all_2[v] <- loss_t_all_2[v] + loss_k
    }
    loss_t_all[v] <- loss_t
    
    if(v>1){
      loss_min <- loss_t_all[v-1] -loss_t_all[v]
    }
  }
  
  loss_final <- loss_t
  stop <- 0
  iter_all <- 0
  while(stop ==0 & iter_all < 50){
    
    iter_all <- iter_all + 1
    member_exchange <- 0
    loss_n_all <- list()
    for(n in 1:all_member){#length(cluster_assign)){
      #loss_n <- rep(NA, n_cluster)
      loss_n_all[[n]] <- rep(NA, n_cluster)
      ## get the cluster membership
      cluster_n <- cluster_assign[n]
      cluster_assign_temp <- cluster_assign
      #data_x_n <- data_x[-n,]
      #data_y_n <- data_y[-n]
      #t_n <- t[-n,]
      #cluster_assign_n <- cluster_assign[-n]
      py_temp <- list()
      px_temp <- list()
      t_temp <- list()
      
      for(g in 1:n_cluster){
        if (g == cluster_n){
          py_temp[[g]] <- py
          px_temp[[g]] <- px
          t_temp[[g]] <- t
          ##################
          loss_n_all[[n]][g] <- loss_final
        }
        if(g != cluster_n){
          cluster_assign_temp[n] <- g
          loss_min <- upper
          py_temp_k <- list()
          #loss_p_all <- rep(NA, iteration)
          loss_t_all <- rep(NA, iteration)
          
          for (v in 1:iteration){
            if(loss_min < converge_2 & v > 10)  break
            loss_new <- rep(NA, all_member)
            
            if (v == 1){
              t_temp_k <- t
              px_temp_k <- t(data_x) %*% t_temp_k
              px_temp_k[distinct_zeros] <- 0
              for(k in 1:n_cluster){
                cluster_k <- which(cluster_assign_temp == k)
                t_k <- t_temp_k[cluster_k,]
                y_k <- data_y[cluster_k,]#-mean(data_y[cluster_k])
                py_temp_k[[k]] <- t(y_k) %*% t_k %*% inv(t(t_k) %*% t_k)
              }
            }
            
            t_temp_k <- matrix(0, nrow = all_member, ncol = n_total)
            for (k in 1:n_cluster){
              cluster_k <- which(cluster_assign_temp == k)
              x_k <- data_x[cluster_k,]
              xp <- beta * (x_k %*% px)
              y_k <- data_y[cluster_k,]# - mean(data_y[cluster_k])
              xp1 <-(1-beta) * (y_k %*% py_temp_k[[k]])
              xp_update <- xp + xp1
              xp_svd <- svd(xp_update)
              t_temp_k[cluster_k,] <- xp_svd$u %*% t(xp_svd$v) / sqrt(n_cluster)
            }
            
            px_temp_k <- t(data_x) %*% t_temp_k
            px_temp_k[distinct_zeros] <- 0
            px_temp[[k]] <- px_temp_k
            loss_t <- beta * sum((data_x - t_temp_k %*% t(px_temp_k))^2)
            for(k in 1:n_cluster){
              cluster_k <- which(cluster_assign_temp == k)
              t_k <- t_temp_k[cluster_k,]
              y_k <- data_y[cluster_k,]# - mean(data_y[cluster_k])
              reg_results <- t(y_k) %*% t_k %*% inv(t(t_k) %*% t_k)
              py_temp_k[[k]]  <- reg_results
              loss_k <- (1-beta) * sum((y_k - t_k %*% t(py_temp_k[[k]]))^2)
              loss_t <- loss_t+loss_k
            }
            loss_t_all[v] <- loss_t
            
            if(v>1){
              loss_min <- loss_t_all[v-1] -loss_t_all[v]
            }
          }
          loss_n_all[[n]][g] <- loss_t
          py_temp[[g]] <- py_temp_k
          px_temp[[g]] <- px_temp_k
          t_temp[[g]] <- t_temp_k
        }
      }
      cluster_assign[n] <- which(loss_n_all[[n]] == min(loss_n_all[[n]]))
      g <- cluster_assign[n]
      py <- py_temp[[g]]
      px <- px_temp[[g]]
      t <- t_temp[[g]]
      loss_final <- min(loss_n_all[[n]])
      
      if(cluster_assign[n] != cluster_n){
        member_exchange <- 1
      }
    }
    
    if (member_exchange == 0){
      stop <- 1
    }
  }
  results <- list(cluster_assign = cluster_assign, loss = loss_n_all, score = t, loadings = px,
                  regs = py)
  return(results)
}

###########################################################################
###########################################################################
### without any model selection
CSSCR_nocon <- function(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster, alpha = .8,
                        converge = 1e-7, converge_2 = 1e-7, iteration = 1000, start_part = NULL){
  
  ## fix the sparsity level to zero (p_sparse = 0)
  upper <- 1e9
  stop <- 0
  iter <- 0
  
  all_member <- nrow(data_x)
  sum_var <- ncol(data_x)
  y_var <- ncol(data_y)
  block_version_data <- list(length = 2)
  block_version_data[[1]] <- data_x[,1:n_var[1]]
  block_version_data[[2]] <- data_x[, (n_var[1] + 1):n_var[2]]
  n_total <- sum(n_com, n_distinct)
  
  ## structure-induced zeros
  distinct_index <- vector("numeric")
  distinct_zeros <- vector("numeric")
  if (sum(ndistinct) != 0){
    
    all_var <- 0
    ini_var <- 0
    for (p in 1:n_block){
      distinct_index <- c(distinct_index, rep(p, n_distinct[p]))
      ini_var <- ini_var + n_var[p]
      all_var <- c(all_var, ini_var)
    }
    
    for (r.distinct in 1:sum(n_distinct)){
      distinct_zeros <- c(distinct_zeros, ((sum_var * (n_com + r.distinct - 1)) + all_var[distinct_index[r.distinct]] + 1): ((sum_var * (n_com + r.distinct - 1)) + all_var[(distinct_index[r.distinct] + 1)]))
    }
  }
  
  # set the upper bound of minimum loss
  loss_all <- vector("numeric", length = iteration)
  
  # the starting partition
  start <- vector("numeric", length = n_cluster)
  for (y in 1:n_cluster){
    start[y] <- round(all_member / n_cluster)
  }
  start[n_cluster] <- all_member - (n_cluster - 1) * round(all_member / n_cluster)
  
  
  py <- list()
  if(is.null(start_part)){
    cluster_assign <- RandomStart(start)[[2]]
  }
  if(!is.null(start_part)){
    cluster_assign <- start_part
  }
  
  ###################################################
  #cluster_assign <- a
  svd_data <- svds(data_x, n_total)
  t <- svd_data$u
  x_square <- sum(data_x^2)
  y_square <- sum(data_y^2)
  #beta <- .5
  #alpha <- .5
  beta <- alpha * y_square / (alpha * y_square + (1-alpha) * x_square)
  
  #data_x <- MatrixCenter(data_x,1,0)
  #beta <- 0.9
  
  loss_min <- upper
  loss_p_all <- rep(NA, iteration)
  loss_t_all <- rep(NA, iteration)
  loss_t_all_1 <- rep(NA, iteration)
  loss_t_all_2 <- rep(NA, iteration)
  loss_t_all_3 <- rep(NA, iteration)
  for (v in 1:iteration){
    if(loss_min < converge & v > 20)  break
    loss_new <- rep(NA, all_member)
    if (v == 1){
      px <- t(data_x) %*% t
      px[distinct_zeros] <- 0
    }
    beta * sum((data_x - t %*% t(px))^2)
    #index.nsmallest <- sort(px^2, index.return = TRUE)$ix[1:(sparse_zeros+length(distinct_zeros))]
    #x[index.nsmallest] <- 0
    
    #if(v > 1){
    #  loss_p <- beta * sum((data_x - t %*% t(px))^2)
    #  for(k in 1:n_cluster){
    #    cluster_k <- which(cluster_assign == k)
    #    t_k <- t[cluster_k,]
    #    y_k <- data_y[cluster_k]
    #    loss_k <- (1-beta) * sum((y_k-(cbind(rep(1, nrow(t_k)), t_k)%*%py[[k]]))^2)
    #    loss_p <- loss_p+loss_k
    #  }
    #  loss_p_all[v] <- loss_p
    #}
    
    if (v == 1){
      for(k in 1:n_cluster){
        cluster_k <- which(cluster_assign == k)
        t_k <- t[cluster_k,]
        y_k <- data_y[cluster_k,]
        py[[k]] <- t(y_k) %*% t_k %*% inv(t(t_k) %*% t_k)
      }
    }
    
    ############################
    
    t <- matrix(0, nrow = all_member, ncol = n_total)
    for (k in 1:n_cluster){
      cluster_k <- which(cluster_assign == k)
      x_k <- data_x[cluster_k, ]
      y_k <- data_y[cluster_k]# - mean(data_y[cluster_k])
      zz <- beta * (diag(length(cluster_k)) %x% (t(px) %*% px)) + (1-beta)*(diag(length(cluster_k)) %x% (t(py[[k]]) %*% py[[k]]))
      zy <- beta * ((diag(length(cluster_k)) %x% t(px)) %*% as.vector(t(x_k))) + (1-beta)*((diag(length(cluster_k)) %x% t(py[[k]])) %*% y_k)
      t[cluster_k, ]<- matrix(inv(zz) %*% zy, nrow = length(cluster_k), byrow = TRUE)
    }
    
    ####################################
    loss_t2 <- 0
    for(k in 1:n_cluster){
      cluster_k <- which(cluster_assign == k)
      t_k <- t[cluster_k,]
      y_k <- data_y[cluster_k,]#-mean(data_y[cluster_k])
      loss_k <- (1-beta) * sum((y_k - t_k %*% t(py[[k]]))^2)
      loss_t2 <- loss_t2+loss_k
    }
    
    px <- t(data_x) %*% t %*% inv(t(t) %*% t)
    px[distinct_zeros] <- 0
    
    ######################################
    loss_t <- beta * sum((data_x - t %*% t(px))^2)
    #loss_t_all_1[v] <- loss_t
    loss_t_all_2[v] <- 0
    for(k in 1:n_cluster){
      cluster_k <- which(cluster_assign == k)
      t_k <- t[cluster_k,]
      y_k <- data_y[cluster_k,]#-mean(data_y[cluster_k])
      reg_results <- t(y_k) %*% t_k %*% inv(t(t_k) %*% t_k)
      py[[k]] <- reg_results
      loss_k <- (1-beta) * sum((y_k - t_k %*% t(py[[k]]))^2)
      loss_t <- loss_t+loss_k
      loss_t_all_2[v] <- loss_t_all_2[v] + loss_k
    }
    loss_t_all[v] <- loss_t
    
    if(v>1){
      loss_min <- loss_t_all[v-1] -loss_t_all[v]
    }
  }
  
  loss_final <- loss_t
  stop <- 0
  iter_all <- 0
  while(stop ==0 & iter_all < 50){
    
    iter_all <- iter_all + 1
    member_exchange <- 0
    loss_n_all <- list()
    for(n in 1:all_member){#length(cluster_assign)){
      #loss_n <- rep(NA, n_cluster)
      loss_n_all[[n]] <- rep(NA, n_cluster)
      ## get the cluster membership
      cluster_n <- cluster_assign[n]
      cluster_assign_temp <- cluster_assign
      #data_x_n <- data_x[-n,]
      #data_y_n <- data_y[-n]
      #t_n <- t[-n,]
      #cluster_assign_n <- cluster_assign[-n]
      py_temp <- list()
      px_temp <- list()
      t_temp <- list()
      
      for(g in 1:n_cluster){
        if (g == cluster_n){
          py_temp[[g]] <- py
          px_temp[[g]] <- px
          t_temp[[g]] <- t
          ##################
          loss_n_all[[n]][g] <- loss_final
        }
        if(g != cluster_n){
          cluster_assign_temp[n] <- g
          loss_min <- upper
          py_temp_k <- list()
          #loss_p_all <- rep(NA, iteration)
          loss_t_all <- rep(NA, iteration)
          
          for (v in 1:iteration){
            if(loss_min < converge_2 & v > 10)  break
            loss_new <- rep(NA, all_member)
            
            if (v == 1){
              t_temp_k <- t
              px_temp_k <- t(data_x) %*% t_temp_k %*% inv(t(t_temp_k) %*% t_temp_k)
              px_temp_k[distinct_zeros] <- 0
              for(k in 1:n_cluster){
                cluster_k <- which(cluster_assign_temp == k)
                t_k <- t_temp_k[cluster_k,]
                y_k <- data_y[cluster_k,]#-mean(data_y[cluster_k])
                py_temp_k[[k]] <- t(y_k) %*% t_k %*% inv(t(t_k) %*% t_k)
              }
            }
            
            
            #loss_p <- beta * sum((data_x - t %*% t(px_temp_k))^2)
            #for(k in 1:n_cluster){
            #  cluster_k <- which(cluster_assign_temp == k)
            #  t_k <- t[cluster_k,]
            #  y_k <- data_y[cluster_k]
            #  loss_k <- (1-beta) * sum((y_k-(cbind(rep(1, nrow(t_k)), t_k)%*%py[[k]]))^2)
            #  loss_p <- loss_p+loss_k
            #}
            #}
            t_temp_k <- matrix(0, nrow = all_member, ncol = n_total)
            for (k in 1:n_cluster){
              cluster_k <- which(cluster_assign_temp == k)
              x_k <- data_x[cluster_k,]
              xp <- beta * (x_k %*% px)
              y_k <- data_y[cluster_k,]# - mean(data_y[cluster_k])
              zz <- beta * (diag(length(cluster_k)) %x% (t(px_temp_k) %*% px_temp_k)) + (1-beta)*(diag(length(cluster_k)) %x% (t(py_temp_k[[k]]) %*% py_temp_k[[k]]))
              zy <- beta * ((diag(length(cluster_k)) %x% t(px_temp_k)) %*% as.vector(t(x_k))) + (1-beta)*((diag(length(cluster_k)) %x% t(py_temp_k[[k]])) %*% y_k)
              t_temp_k[cluster_k,] <- matrix(inv(zz) %*% zy, nrow = length(cluster_k), byrow = TRUE)
            }
            
            
            #t.ee <- eigen(t(t_new) %*% t_new, symmetric = FALSE)
            #t <- t_new %*% t.ee$vectors %*% inv(sqrt(diag(t.ee$values)))
            #gs <- gramSchmidt(t_new)
            #q <- px_temp_k %*% t(gs$R)
            #t_temp_k <- gs$Q
            
            px_temp_k <- t(data_x) %*% t_temp_k %*% inv(t(t_temp_k) %*% t_temp_k)
            px_temp_k[distinct_zeros] <- 0
            px_temp[[k]] <- px_temp_k
            loss_t <- beta * sum((data_x - t_temp_k %*% t(px_temp_k))^2)
            for(k in 1:n_cluster){
              cluster_k <- which(cluster_assign_temp == k)
              t_k <- t_temp_k[cluster_k,]
              y_k <- data_y[cluster_k,]# - mean(data_y[cluster_k])
              reg_results <- t(y_k) %*% t_k %*% inv(t(t_k) %*% t_k)
              py_temp_k[[k]]  <- reg_results
              loss_k <- (1-beta) * sum((y_k - t_k %*% t(py_temp_k[[k]]))^2)
              loss_t <- loss_t+loss_k
            }
            loss_t_all[v] <- loss_t
            
            if(v>1){
              loss_min <- loss_t_all[v-1] -loss_t_all[v]
            }
          }
          loss_n_all[[n]][g] <- loss_t
          py_temp[[g]] <- py_temp_k
          px_temp[[g]] <- px_temp_k
          t_temp[[g]] <- t_temp_k
        }
      }
      cluster_assign[n] <- which(loss_n_all[[n]] == min(loss_n_all[[n]]))
      g <- cluster_assign[n]
      py <- py_temp[[g]]
      px <- px_temp[[g]]
      t <- t_temp[[g]]
      loss_final <- min(loss_n_all[[n]])
      
      if(cluster_assign[n] != cluster_n){
        member_exchange <- 1
      }
    }
    
    if (member_exchange == 0){
      stop <- 1
    }
  }
  results <- list(cluster_assign = cluster_assign, loss = loss_n_all, score = t, loadings = px,
                  regs = py)
  return(results)
}
#########################################################################
#########################################################################
#########################################################################
MultiCSCCR <- function(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster, alpha = .8,
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


###################
### without any model selection
CSSCR_group <- function(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster,  alpha = .8,
                  converge = 1e-8, converge_2 = 1e-8, iteration = 1000, start_part = NULL){
  
  ## fix the sparsity level to zero (p_sparse = 0)
  upper <- 1e9
  stop <- 0
  iter <- 0
  
  all_member <- nrow(data_x)
  sum_var <- ncol(data_x)
  y_var <- ncol(data_y)
  block_version_data <- list(length = 2)
  block_version_data[[1]] <- data_x[,1:n_var[1]]
  block_version_data[[2]] <- data_x[, (n_var[1] + 1):n_var[2]]
  n_total <- sum(n_com, n_distinct)

  
  ## structure-induced zeros
  distinct_index <- vector("numeric")
  distinct_zeros <- vector("numeric")
  if (sum(n_distinct) != 0){
    
    all_var <- 0
    ini_var <- 0
    for (p in 1:n_block){
      distinct_index <- c(distinct_index, rep(p, n_distinct[p]))
      ini_var <- ini_var + n_var[p]
      all_var <- c(all_var, ini_var)
    }
    
    for (r.distinct in 1:sum(n_distinct)){
      distinct_zeros <- c(distinct_zeros, ((sum_var * (n_com + r.distinct - 1)) + all_var[distinct_index[r.distinct]] + 1): ((sum_var * (n_com + r.distinct - 1)) + all_var[(distinct_index[r.distinct] + 1)]))
    }
  }
  
  # set the upper bound of minimum loss
  loss_all <- vector("numeric", length = iteration)
  
  # the starting partition
  start <- vector("numeric", length = n_cluster)
  for (y in 1:n_cluster){
    start[y] <- round(all_member / n_cluster)
  }
  start[n_cluster] <- all_member - (n_cluster - 1) * round(all_member / n_cluster)
  
  
  py <- list()
  px <- list()
  
  if(is.null(start_part)){
    cluster_assign <- RandomStart(start)[[2]]
  }
  if(!is.null(start_part)){
    cluster_assign <- start_part
  }
  
  ###################################################
  #cluster_assign <- a
  data_x <- MatrixCenter(data_x,1,0)
  data_y <- data_y - mean(data_y)
  svd_data <- svd(data_x, n_total)
  t <- svd_data$u
  x_square <- sum(data_x^2)
  y_square <- sum(data_y^2)
  #beta <- .5
  #alpha <- .5
  beta <- alpha * y_square / (alpha * y_square + (1-alpha) * x_square)
  
  #data_x <- MatrixCenter(data_x,1,0)
  #beta <- 0.9
  
  loss_min <- upper
  loss_p_all <- rep(NA, iteration)
  loss_t_all <- rep(NA, iteration)
  #loss_t_all_1 <- rep(NA, iteration)
  #loss_t_all_2 <- rep(NA, iteration)
  #loss_t_all_3 <- rep(NA, iteration)
  mean_x <- matrix(nrow = n_cluster, ncol = sum_var)
  mean_y <- rep(NA, n_cluster)
  for (v in 1:iteration){
    if(loss_min < converge & v > 20)  break
    loss_new <- rep(NA, all_member)
    
    if (v == 1){
      for(k in 1:n_cluster){
        cluster_k <- which(cluster_assign == k)
        t_k <- t[cluster_k,]
        x_k <- data_x[cluster_k,]
        x_k <- MatrixCenter(x_k,1,0)
        px[[k]] <- t(x_k) %*% t_k
        px[[k]][distinct_zeros] <- 0
        y_k <- data_y[cluster_k]# - mean(data_y[cluster_k])
        t_k_comb <- cbind(rep(1,nrow(t_k)), t_k)
        py[[k]] <- t(y_k) %*% t_k_comb %*% inv(t(t_k_comb) %*% t_k_comb)
      }
    }
    
    ############################
    
    t <- matrix(0, nrow = all_member, ncol = n_total)
    for (k in 1:n_cluster){
      cluster_k <- which(cluster_assign == k)
      x_k <- data_x[cluster_k,]
      x_k <- MatrixCenter(x_k,1,0)
      xp <- beta * (x_k %*% px[[k]])
      y_k <- data_y[cluster_k] - py[[k]][1]# - mean(data_y[cluster_k])
      xp1 <- (1-beta) * (y_k %*% t(py[[k]][-1]))
      xp_update <- xp + xp1
      xp_svd <- svd(xp_update)
      t[cluster_k,] <- (xp_svd$u %*% t(xp_svd$v)) 
    }
    ####################################
    #px <- t(data_x) %*% t
    #px[distinct_zeros] <- 0
    
    #loss_t_all_1[v] <- loss_t
    loss_t_all[v] <- 0
    for(k in 1:n_cluster){
      cluster_k <- which(cluster_assign == k)
      t_k <- t[cluster_k,]
      x_k <- data_x[cluster_k,]
      mean_x[k,] <- apply(x_k, 2, mean)
      x_k <- MatrixCenter(x_k,1,0)
      px[[k]] <- t(x_k) %*% t_k
      px[[k]][distinct_zeros] <- 0
      y_k <- data_y[cluster_k]#-mean(data_y[cluster_k])
      #mean_y[k] <- mean(data_y[cluster_k])
      t_k_comb <- cbind(rep(1,nrow(t_k)), t_k)
      reg_results <- t(y_k) %*% t_k_comb %*% inv(t(t_k_comb) %*% t_k_comb)
      py[[k]] <- reg_results
      loss_k <- (1-beta) * sum((y_k - t_k_comb %*% t(py[[k]]))^2) + beta * sum((x_k - t_k %*% t(px[[k]]))^2)
      loss_t_all[v] <- loss_t_all[v] + loss_k
    }
    
    if(v>1){
      loss_min <- loss_t_all[v-1] -loss_t_all[v]
    }
  }
  
  loss_final <- loss_t_all[v-1]
  stop <- 0
  iter_all <- 0
  while(stop ==0 & iter_all < 50){
    
    iter_all <- iter_all + 1
    member_exchange <- 0
    loss_n_all <- list()
    for(n in 1:all_member){#length(cluster_assign)){
      #loss_n <- rep(NA, n_cluster)
      loss_n_all[[n]] <- rep(NA, n_cluster)
      ## get the cluster membership
      cluster_n <- cluster_assign[n]
      cluster_assign_temp <- cluster_assign
      #data_x_n <- data_x[-n,]
      #data_y_n <- data_y[-n]
      #t_n <- t[-n,]
      #cluster_assign_n <- cluster_assign[-n]
      py_temp <- list()
      px_temp <- list()
      t_temp <- list()
      mean_x_temp <- list()
      mean_y_temp <- list()
      
      for(g in 1:n_cluster){
        if (g == cluster_n){
          py_temp[[g]] <- py
          px_temp[[g]] <- px
          t_temp[[g]] <- t
          mean_x_temp[[g]] <- mean_x
          mean_y_temp[[g]] <- mean_y
          ##################
          loss_n_all[[n]][g] <- loss_final
        }
        if(g != cluster_n){
          cluster_assign_temp[n] <- g
          loss_min <- upper
          py_temp_k <- list()
          px_temp_k <- list()
          #loss_p_all <- rep(NA, iteration)
          loss_t_all <- rep(NA, iteration)
          mean_x_temp[[g]] <- mean_x
          mean_y_temp[[g]] <- mean_y
          
          for (v in 1:iteration){
            if(loss_min < converge_2 & v > 10)  break
            loss_new <- rep(NA, all_member)
            
            if (v == 1){
              t_temp_k <- t
              for(k in 1:n_cluster){
                cluster_k <- which(cluster_assign_temp == k)
                x_k <- data_x[cluster_k,]
                x_k <- MatrixCenter(x_k,1,0)
                t_k <- t_temp_k[cluster_k,]
                px_temp_k[[k]] <- t(x_k) %*% t_k
                px_temp_k[[k]][distinct_zeros] <- 0
                y_k <- data_y[cluster_k]#-mean(data_y[cluster_k])
                t_k_comb <- cbind(rep(1,nrow(t_k)), t_k)
                py_temp_k[[k]] <- t(y_k) %*% t_k_comb %*% inv(t(t_k_comb) %*% t_k_comb)
              }
            }
            
            t_temp_k <- matrix(0, nrow = all_member, ncol = n_total)
            for (k in 1:n_cluster){
              cluster_k <- which(cluster_assign_temp == k)
              x_k <- data_x[cluster_k,]

              x_k <- MatrixCenter(x_k,1,0)
              xp <- beta * (x_k %*% px_temp_k[[k]])
              y_k <- data_y[cluster_k]-py_temp_k[[k]][1]# - mean(data_y[cluster_k])

              xp1 <-(1-beta) * (y_k %*% t(py_temp_k[[k]][-1]))
              xp_update <- xp + xp1
              xp_svd <- svd(xp_update)
              t_temp_k[cluster_k,] <- xp_svd$u %*% t(xp_svd$v)# / sqrt(n_cluster)
            }
            
            loss_t <- 0
            for(k in 1:n_cluster){
              cluster_k <- which(cluster_assign_temp == k)
              t_k <- t_temp_k[cluster_k,]
              y_k <- data_y[cluster_k]# - mean(data_y[cluster_k])
              mean_y_temp[[g]][k] <- mean(data_y[cluster_k])
              x_k <- data_x[cluster_k,]
              mean_x_temp[[g]][k,] <- apply(x_k,2,mean)
              x_k <- MatrixCenter(x_k,1,0)
              px_temp_k[[k]] <- t(x_k) %*% t_k
              px_temp_k[[k]][distinct_zeros] <- 0
              t_k_comb <- cbind(rep(1,nrow(t_k)), t_k) 
              reg_results <- t(y_k) %*% t_k_comb %*% inv(t(t_k_comb) %*% t_k_comb)
              py_temp_k[[k]]  <- reg_results
              loss_k <- (1-beta) * sum((y_k - t_k_comb %*% t(py_temp_k[[k]]))^2) + beta * sum((x_k - t_k %*% t(px_temp_k[[k]]))^2)
              loss_t <- loss_t+loss_k
            }
            loss_t_all[v] <- loss_t
            
            if(v>1){
              loss_min <- loss_t_all[v-1] -loss_t_all[v]
            }
          }
          loss_n_all[[n]][g] <- loss_t
          py_temp[[g]] <- py_temp_k
          px_temp[[g]] <- px_temp_k
          t_temp[[g]] <- t_temp_k
        }
      }
      cluster_assign[n] <- which(loss_n_all[[n]] == min(loss_n_all[[n]]))
      g <- cluster_assign[n]
      py <- py_temp[[g]]
      px <- px_temp[[g]]
      t <- t_temp[[g]]
      mean_x <- mean_x_temp[[g]] 
      mean_y <- mean_y_temp[[g]]
      loss_final <- min(loss_n_all[[n]])
      
      if(cluster_assign[n] != cluster_n){
        member_exchange <- 1
      }
    }
    
    if (member_exchange == 0){
      stop <- 1
    }
  }
  results <- list(cluster_assign = cluster_assign, loss = loss_final, score = t, loadings = px,
                  regs = py, loss_all = loss_n_all, loss_min = loss_min, mean_x = mean_x, mean_y = mean_y)
  return(results)
}

CSSCRPredict <- function(old_data_x=NULL, new_data_x=NULL, old_loading_x, old_loading_y, mean_x, new_assign = NULL){
  
  n_new <- nrow(new_data_x)
  n_cluster <- length(old_loading_x)
  n_var <- nrow(old_loading_x[[1]])
  
  if(is.null(new_assign)){
    loss_x <- matrix(nrow = n_new, ncol = n_cluster)
    for (i in 1:n_new){
      for(j in 1:n_cluster){
        score_i <- t(new_data_x[i,]- mean_x[j,]) %*% old_loading_x[[j]] %*% inv(t(old_loading_x[[j]]) %*% old_loading_x[[j]])
        loss_x[i,j] <- sum((new_data_x[i,] - mean_x[j,]- score_i %*% t(old_loading_x[[j]]))^2)
      }
    }
    new_assign <- apply(loss_x,1,function(x) which.min(x))
  }
  predict_y <- rep(NA,n_new) 
  for(i in 1:n_new){
    new_cluster <- new_assign[i]
    loading_i <- old_loading_y[[new_cluster]]
    score_i <- t(new_data_x[i,]- mean_x[new_cluster,]) %*% old_loading_x[[new_cluster]] %*% inv(t(old_loading_x[[new_cluster]]) %*% old_loading_x[[new_cluster]])
    score_i_new <- cbind(rep(1,nrow(score_i)), score_i)
    predict_y[i] <- score_i_new %*% t(loading_i)
  }
  
  results <- list(predict_y = predict_y, predict_cluster = new_assign)
  return(results)
}

#######################################
###################
### without any model selection
CSCCR_single <- function(data_x_complete,data_y_complete,n_block, n_com, n_distinct, n_var, alpha = .8,
                        converge = 1e-8, converge_2 = 1e-8, iteration = 1000, cluster_assign){
  
  ## fix the sparsity level to zero (p_sparse = 0)
  upper <- 1e9
  stop <- 0
  iter <- 0
  
  all_member <- nrow(data_x)
  sum_var <- ncol(data_x)
  y_var <- ncol(data_y)
  block_version_data <- list(length = 2)
  block_version_data[[1]] <- data_x[,1:n_var[1]]
  block_version_data[[2]] <- data_x[, (n_var[1] + 1):n_var[2]]
  n_total <- sum(n_com, n_distinct)
  
  
  ## structure-induced zeros
  distinct_index <- vector("numeric")
  distinct_zeros <- vector("numeric")
  if (sum(ndistinct) != 0){
    
    all_var <- 0
    ini_var <- 0
    for (p in 1:n_block){
      distinct_index <- c(distinct_index, rep(p, n_distinct[p]))
      ini_var <- ini_var + n_var[p]
      all_var <- c(all_var, ini_var)
    }
    
    for (r.distinct in 1:sum(n_distinct)){
      distinct_zeros <- c(distinct_zeros, ((sum_var * (n_com + r.distinct - 1)) + all_var[distinct_index[r.distinct]] + 1): ((sum_var * (n_com + r.distinct - 1)) + all_var[(distinct_index[r.distinct] + 1)]))
    }
  }
  
  # set the upper bound of minimum loss
  loss_all <- vector("numeric", length = iteration)
  ###################################################
  px <- list()
  py <- list()
  t <- matrix(nrow = all_member, ncol = n_total)
  #cluster_assign <- a
  n_cluster <- length(unique(cluster_assign))
  mean_x <- matrix(nrow = n_cluster, ncol = sum_var)
  for (i in 1:n_cluster){
    data_x <- data_x_complete[which(cluster_assign == i),]
    mean_x[i,] <- apply(data_x,2,mean)
    data_x <- MatrixCenter(data_x,1,0)

    data_y <- data_y_complete[which(cluster_assign == i)]
    data_y_center <- data_y - mean(data_y)
    svd_data <- svd(data_x, n_total)
    t_i <- svd_data$u
    x_square <- sum(data_x^2)
    y_square <- sum(data_y_center^2)
    #beta <- .5
    #alpha <- .5
    beta <- alpha * y_square / (alpha * y_square + (1-alpha) * x_square)
    
    #data_x <- MatrixCenter(data_x,1,0)
    #beta <- 0.9
    
    loss_min <- upper
    loss_t_all <- rep(NA, iteration)
  
    for (v in 1:iteration){
      if(loss_min < converge & v > 20)  break
      loss_new <- rep(NA, all_member)
      
      if (v == 1){
  
          px_i <- t(data_x) %*% t_i
          px_i[distinct_zeros] <- 0
          t_comb <- cbind(rep(1,nrow(t_i)), t_i)
          py_i <- t(data_y) %*% t_comb %*% inv(t(t_comb) %*% t_comb)
      }
      
      ############################
  
        xp <- beta * (data_x %*% px_i)
        y_k <- data_y - py_i[1]# - mean(data_y[cluster_k])
        xp1 <- (1-beta) * (y_k %*% t(py_i[-1]))
        xp_update <- xp + xp1
        xp_svd <- svd(xp_update)
        t_i <- (xp_svd$u %*% t(xp_svd$v)) 
      ####################################
      #px <- t(data_x) %*% t
      #px[distinct_zeros] <- 0
      
      #loss_t_all_1[v] <- loss_t
      loss_t_all[v] <- 0
  
      px_i <- t(data_x) %*% t_i
      px_i[distinct_zeros] <- 0
      t_comb <- cbind(rep(1,nrow(t_i)), t_i)
      py_i <- t(data_y) %*% t_comb %*% inv(t(t_comb) %*% t_comb)
      loss_k <- (1-beta) * sum((data_y - t_comb %*% t(py_i))^2) + beta * sum((data_x - t_i %*% t(px_i))^2)
      loss_t_all[v] <- loss_t_all[v] + loss_k
      
      if(v>1){
        loss_min <- loss_t_all[v-1] -loss_t_all[v]
      }
    }
    
    t[which(cluster_assign == i),] <- t_i
    px[[i]] <- px_i
    py[[i]] <- py_i
  }
  
  
  results <- list(loss = loss_min, score = t, loadings = px,
                  regs = py, mean_x = mean_x)
  return(results)
}
###############################################################
###############################################################
## the following code is not used
###############################################################


















####################################################################################
####################################################################################
### without any model selection
CSSCRnon <- function(data_x,data_y,n_block, n_com, n_distinct, n_var, n_cluster, p_sparse = 0, beta,
                     converge = 1e-9, iteration = 1000){

  ## fix the sparsity level to zero (p_sparse = 0)
  upper <- 1e9
  stop <- 0
  iter <- 0
  converge_2 <- 1e-6
  converge <- 1e6

  all_member <- nrow(data_x)
  sum_var <- ncol(data_x)
  block_version_data <- list(length = 2)
  block_version_data[[1]] <- data_x[,1:n_var[1]]
  block_version_data[[2]] <- data_x[, (n_var[1] + 1):n_var[2]]
  n_total <- sum(n_com, n_distinct)

  ## structure-induced zeros
  distinct_index <- vector("numeric")
  all_var <- 0
  ini_var <- 0
  for (p in 1:n_block){
    distinct_index <- c(distinct_index, rep(p, n_distinct[p]))
    ini_var <- ini_var + n_var[p]
    all_var <- c(all_var, ini_var)
  }
  distinct_zeros <- vector("numeric")
  for (r.distinct in 1:sum(n_distinct)){
    distinct_zeros <- c(distinct_zeros, ((sum_var * (n_com + r.distinct - 1)) + all_var[distinct_index[r.distinct]] + 1): ((sum_var * (n_com + r.distinct - 1)) + all_var[(distinct_index[r.distinct] + 1)]))
  }

  ## sparsity-induced zeros
  # the positions that are not yet specified as zeros
  #retain_zeros <- setdiff(1: (sum_var * n_total), distinct_zeros)

  # generate the component loading matrix for the predictors
  # the number of zeros in component loadings
  #sparse_zeros <- round(length(retain_zeros) * psparse)

  # set the upper bound of minimum loss
  loss_all <- vector("numeric", length = iteration)

  # the starting partition
  start <- vector("numeric", length = n_cluster)
  for (y in 1:n_cluster){
    start[y] <- round(all_member / n_cluster)
  }
  start[n_cluster] <- all_member - (n_cluster - 1) * round(all_member / n_cluster)

  py <- list()
  cluster_assign <- RandomStart(start)[[2]]
  ###################################################
  cluster_assign <- c(rep(1,25), rep(2,25), rep(3,25))
  cluster_assign[1:5]<- 3
  cluster_assign[71:75] <- 1
  cluster_assign <- rep(c(rep(1,10),rep(2,10), rep(3,10)),5)
  ###################################################
  #cluster_assign <- a
  t_matrix <- matrix(nrow = all_member, ncol = n_total)
  x_square <- sum(data_x^2)
  y_square <- sum(data_y^2)
  alpha <- .1
  beta <- alpha * y_square / (alpha * y_square + (1-alpha) * x_square)
  #beta <- 0.1

  loss_p_all <- rep(NA, iteration)
  loss_t_all <- 0
  loss_t_all_1 <- list()
  px <- list()

  for(k in 1:n_cluster){
    loss_t_all_1[[k]] <- rep(NA, iteration)
    cluster_k <- which(cluster_assign == k)
    data_x_k <- data_x[cluster_k, ]
    data_x_k_c <- MatrixCenter(data_x_k, 1, 0)
    data_y_k <- data_y[cluster_k]
    data_y_k_c <- data_y[cluster_k] - mean(data_y[cluster_k])
    svd_data <- svds(data_x_k_c, n_total)
    t_k <- svd_data$u
    loss_min <- upper

    for (v in 1:iteration){
      if(loss_min < converge & v > 20)  break
      if (v == 1){
        px[[k]] <- t(data_x_k_c) %*% t_k
        px[[k]][distinct_zeros] <- 0
      }
      reg_results <- lm(data_y_k_c~t_k-1)
      py[[k]] <- reg_results$coefficients
      t_new <- t_k
      for(i in 1:length(cluster_k)){
        #update_i <- which(cluster_k == i)
        q <- py[[k]]
        yi <- data_y_k_c[i]
        ti <- t(inv(beta * q %*% t(q) - q %*% t(q) - beta * t(px[[k]]) %*% px[[k]]) %*% t((beta*yi*t(q) - yi*t(q) - beta * data_x_k_c[i,] %*% px[[k]])))
        t_new[i,] <- ti
      }

      t.de <- svd(t_new)
      t_new_ortho <-  t.de$u %*% t(t.de$v)
      t_k <- t_new_ortho

      px[[k]] <- t(data_x_k_c) %*% t_k
    #pp <- px
      px[[k]][distinct_zeros] <- 0
      loss_t <- beta * sum((data_x_k_c - t_k %*% t(px[[k]]))^2)
      loss_t_all_1[[k]][v] <- loss_t

      if(v>1){
        loss_min <- loss_t_all_1[[k]][v-1] -loss_t_all_1[[k]][v]
      }
    }
    t_matrix[cluster_k,] <- t_k
    loss_t_all <- loss_t_all + loss_t
  }

  #loss_t_all_1[[3]][1:10]
    #for(k in 1:n_cluster){
    #  cluster_k <- which(cluster_assign == k)
    #  t_k <- t[cluster_k,]
    #  y_k <- data_y[cluster_k]-mean(data_y[cluster_k])
    ##  reg_results <- lm(y_k~t_k-1)
    #  py[[k]] <- reg_results$coefficients
    #  loss_k <- (1-beta) * sum(reg_results$residuals^2)
    #  loss_t <- loss_t+loss_k
    #  loss_t_all_1[v] <- loss_t_all_1[v] + loss_k
    #}
    #loss_t_all[v] <- loss_t
  #}

  #stop <- 0
  #iter_all <- 0
  #while(stop ==0){

   # iter_all <- iter_all + 1
    member_exchange <- 0
    loss_n_all <- list()
    for(n in 1:75){#:10){#length(cluster_assign)){
      #loss_n <- rep(NA, n_cluster)
      loss_n_all[[n]] <- rep(NA, n_cluster)
      ## get the cluster membership
      cluster_n <- cluster_assign[n]
      cluster_assign_temp <- cluster_assign
      #data_x_n <- data_x[-n,]
      #data_y_n <- data_y[-n]
      #t_n <- t[-n,]
      #cluster_assign_n <- cluster_assign[-n]
      py_temp <- list()
      px_temp <- list()
      t_temp <- list()

      for(g in 1:n_cluster){
        if (g == cluster_n){
          py_temp[[g]] <- py
          px_temp[[g]] <- px
          t_temp[[g]] <- t_matrix
          ##################
          #px <- t(data_x) %*% t
          #px[distinct_zeros] <- 0
          #loss_t <- beta * sum((data_x - t %*% t(px))^2)
          #for(k in 1:n_cluster){
          #  cluster_k <- which(cluster_assign == k)
          #  t_k <- t[cluster_k,]
          #
           # y_k <- data_y[cluster_k] - mean(data_y[cluster_k])

            #reg_results <- lm(y_k~t_k-1)
            #py[[k]] <- reg_results$coefficients
            #loss_k <- (1-beta) * sum(reg_results$residuals^2)
            #loss_t <- loss_t+loss_k
          #}
          loss_n_all[[n]][g] <- loss_t_all
        }
        if(g != cluster_n){
          cluster_assign_temp[n] <- g
          py_temp[[g]] <- list()
          px_temp[[g]] <- list()
          t_temp[[g]] <- matrix(nrow = nrow(t_matrix), ncol = ncol(t_matrix))
          loss_n_all[[n]][g] <- 0
          for(k in 1:n_cluster){
            loss_t_all_1[[k]] <- rep(NA, iteration)
            cluster_k <- which(cluster_assign_temp == k)
            data_x_k <- data_x[cluster_k, ]
            data_x_k_c <- MatrixCenter(data_x_k, 1, 0)
            data_y_k <- data_y[cluster_k]
            data_y_k_c <- data_y[cluster_k] - mean(data_y[cluster_k])
            t_k <- t_matrix[cluster_k,]
            loss_min <- upper

            for (v in 1:iteration){
              if(loss_min < converge & v > 20)  break
              if (v == 1){
                px_temp[[g]][[k]] <- t(data_x_k_c) %*% t_k
                px_temp[[g]][[k]][distinct_zeros] <- 0
              }
              reg_results <- lm(data_y_k_c~t_k-1)
              py_temp[[g]][[k]] <- reg_results$coefficients
              t_new <- t_k
              for(i in 1:length(cluster_k)){
                #update_i <- which(cluster_k == i)
                q <- py_temp[[g]][[k]]
                yi <- data_y_k_c[i]
                ti <- t(inv(beta * q %*% t(q) - q %*% t(q) - beta * t(px_temp[[g]][[k]]) %*% px_temp[[g]][[k]]) %*% t((beta*yi*t(q) - yi*t(q) - beta * data_x_k_c[i,] %*% px_temp[[g]][[k]])))
                t_new[i,] <- ti
              }

              t.de <- svd(t_new)
              t_new_ortho <-  t.de$u %*% t(t.de$v)
              t_k <- t_new_ortho

              px_temp[[g]][[k]] <- t(data_x_k_c) %*% t_k
              #pp <- px
              px_temp[[g]][[k]][distinct_zeros] <- 0
              loss_t <- beta * sum((data_x_k_c - t_k %*% t(px_temp[[g]][[k]]))^2)
              loss_t_all_1[[k]][v] <- loss_t

              if(v>1){
                loss_min <- loss_t_all_1[[k]][v-1] -loss_t_all_1[[k]][v]
              }
            }
            t_temp[[g]][cluster_k,] <- t_k
            loss_n_all[[n]][g] <- loss_n_all[[n]][g] + loss_t
          }
        }
      }

        cluster_assign[n] <- which(loss_n_all[[n]] == min(loss_n_all[[n]]))
        g <- cluster_assign[n]
        py <- py_temp[[g]]
        px <- px_temp[[g]]
        t_matrix <- t_temp[[g]]
        loss_t_all <- min(loss_n_all[[n]])
        if(cluster_assign[n] != cluster_n){
          member_exchange <- 1
        }
    }
    #if (member_exchange == 0){
    #  stop <- 1
    #}
    #if (iter_all > 100){
    #  stop <- 1
    #}
  #}





########################################
  ########################
      if (cluster_assign[n] != cluster_n){
        member_exchange <- 1
        #### the update of all elements
        loss_min <- upper
        loss_p_all <- rep(NA, iteration)
        loss_t_all <- rep(NA, iteration)
        for (v in 1:iteration){
          if(loss_min < converge)  break
          loss_new <- rep(NA, all_member)
          px <- t(data_x) %*% t
          px[distinct_zeros] <- 0
          #index.nsmallest <- sort(px^2, index.return = TRUE)$ix[1:(sparse_zeros+length(distinct_zeros))]
          #px[index.nsmallest] <- 0
          if(v > 1){
            loss_p <- beta * sum((data_x - t %*% t(px))^2)
            for(k in 1:n_cluster){
              cluster_k <- which(cluster_assign == k)
              t_k <- t[cluster_k,]
              y_k <- data_y[cluster_k]
              loss_k <- (1-beta) * sum((y_k-(cbind(rep(1, nrow(t_k)), t_k)%*%py[[k]]))^2)
              loss_p <- loss_p+loss_k
            }
            loss_p_all[v] <- loss_p
          }

          t_new <- t
          for(k in 1:n_cluster){
            cluster_k <- which(cluster_assign == k)
            t_k <- t[cluster_k,]
            y_k <- data_y[cluster_k]
            reg_results <- lm(y_k~t_k)
            py[[k]] <- reg_results$coefficients
            for(i in 1:length(cluster_k)){
              update_i <- cluster_k[i]
              q <- py[[k]][-1]
              yi <- data_y[update_i] - py[[k]][1]
              ti <- t(inv(beta * q %*% t(q) - q %*% t(q) - beta * t(px) %*% px) %*% t((beta*yi*t(q) - yi*t(q) - beta * data_x[update_i,] %*% px)))
              t_new[update_i,] <- ti
            }
          }

          t.de <- gramSchmidt(t_new)
          t <- t.de$Q

          px <- t(data_x) %*% t
          px[distinct_zeros] <- 0
          loss_t <- beta * sum((data_x - t %*% t(px))^2)
          for(k in 1:n_cluster){
            cluster_k <- which(cluster_assign == k)
            t_k <- t[cluster_k,]
            y_k <- data_y[cluster_k]
            reg_results <- lm(y_k~t_k)
            py[[k]] <- reg_results$coefficients
            loss_k <- (1-beta) * sum((y_k-(cbind(rep(1, nrow(t_k)), t_k)%*%py[[k]]))^2)
            loss_t <- loss_t+loss_k
          }
          loss_t_all[v] <- loss_t

          if(v>1){
            loss_min <- loss_t_all[v-1] -loss_t_all[v]
          }
        }
      }
    }


  #cluster_assign <- cluster_assign_new

  #}




    member_exchange <- 0
    iter <- iter + 1




      #if(n == 1){
      #  t_new_all <- rbind(t[n,], t_new)
      #}
      #if(n == n_observation){
      #  t_new_all <- rbind(t_new, t[n,])
      #}
      #if(n != 1 & n != n_observation){
      #  t_new_all <- rbind(t_new[1:(n-1),], t[n,], t_new[(n+1):nrow(t_new),])
      #}

      #t.de <- svd(t_new_all)
      #t_new_ortho <-  t.de$u %*% t(t.de$v)
      #t <- t_new_ortho


    if (member_exchange == 0){
      stop <- 1
    }
    if (iter == iteration){
      stop <- 1
    }
    #### the final update
    loss_min <- upper
    loss_p_all <- rep(NA, iteration)
    loss_t_all <- rep(NA, iteration)
    for (v in 1:iteration){
      if(loss_min < converge)  break
      loss_new <- rep(NA, all_member)
      px <- t(data_x) %*% t
      px[distinct_zeros] <- 0
      index.nsmallest <- sort(px^2, index.return = TRUE)$ix[1:(sparse_zeros+length(distinct_zeros))]
      px[index.nsmallest] <- 0
      if(v > 1){
        loss_p <- beta * sum((data_x - t %*% t(px))^2)
        for(k in 1:n_cluster){
          cluster_k <- which(cluster_assign == k)
          t_k <- t[cluster_k,]
          y_k <- data_y[cluster_k]
          loss_k <- (1-beta) * sum((y_k-(cbind(rep(1, nrow(t_k)), t_k)%*%py[[k]]))^2)
          loss_p <- loss_p+loss_k
        }
        loss_p_all[v] <- loss_p
      }

      t_new <- t
      for(k in 1:n_cluster){
        cluster_k <- which(cluster_assign == k)
        t_k <- t[cluster_k,]
        y_k <- data_y[cluster_k]
        reg_results <- lm(y_k~t_k)
        py[[k]] <- reg_results$coefficients
        for(i in 1:length(cluster_k)){
          update_i <- cluster_k[i]
          q <- py[[k]][-1]
          yi <- data_y[update_i] - py[[k]][1]
          ti <- t(inv(beta * q %*% t(q) - q %*% t(q) - beta * t(px) %*% px) %*% t((beta*yi*t(q) - yi*t(q) - beta * data_x[update_i,] %*% px)))
          t_new[update_i,] <- ti
        }
      }

      t.de <- svd(t_new)
      t_new_ortho <-  t.de$u %*% t(t.de$v)
      t <- t_new_ortho

      loss_t <- beta * sum((data_x - t %*% t(px))^2)
      for(k in 1:n_cluster){
        cluster_k <- which(cluster_assign == k)
        t_k <- t[cluster_k,]
        y_k <- data_y[cluster_k]
        reg_results <- lm(y_k~t_k)
        py[[k]] <- reg_results$coefficients
        loss_k <- (1-beta) * sum((y_k-(cbind(rep(1, nrow(t_k)), t_k)%*%py[[k]]))^2)
        loss_t <- loss_t+loss_k
      }
      loss_t_all[v] <- loss_t

      if(v>1){
        loss_min <- loss_t_all[v-1] -loss_t_all[v]
      }
    }
  }
  results = list(cluster_assign = cluster_assign, loadings = px, reg_coef = py, loss = loss_t_all[length(loss_t_all)])
  return(results)
}
