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

  svd_data <- svds(data_x, n_total)
  t <- svd_data$u
  x_square <- sum(data_x^2)
  y_square <- sum(data_y^2)
  beta <- alpha * y_square / (alpha * y_square + (1-alpha) * x_square)

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
