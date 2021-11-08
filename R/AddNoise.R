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
