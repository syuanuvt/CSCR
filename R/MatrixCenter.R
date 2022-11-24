#' Data pre-processing
#'
#' Center and standardized each column of a data matrix
#'
#' @param matrix the matrix to be processed
#' @param center whether column-wise centering is to be carried out
#' @param scale whether column-wise standardization is to be carried out
#' @return A processed matrix
#' @export
#' 
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