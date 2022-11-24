#' Make prediction with the CSCR model
#' 
#' This function allows modellers to predict an new observations based on the information with the existing CSSCR analyis
#'
#' @param old_data_x a matrix indicating the training predictor set
#' @param new_data_x a vector indicating the test predictor set of a new observation
#' @param old_loading_x a matrix indicating the loading matrix associated with the training data set
#' @param old_loading_y  a vector indicating the regression coefficients associated with the training data set
#' @param mean_x the proportion of cluster differences explained by the mean-level differences
#' @param new_assign whether the cluster assignment of the new observation is known
#' @return two items indicate the prediction of the outcome and of the cluster assignment for the new entry
#' @export
#'
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
