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
