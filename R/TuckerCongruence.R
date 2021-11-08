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
