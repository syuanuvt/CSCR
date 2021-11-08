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
