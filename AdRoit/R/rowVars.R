#' Compute row-wise variance of a matrix or data.frame
#'
#' The function computes variance for each row of a matrix or data.frame.
#'
#' @param X a matrix or data.frame.
#' @return a vector of row variances
#' @export

rowVars <- function (X){
  sqr = function(y) y * y
  n = rowSums(!is.na(X))
  n[n <= 1] = NA
  return(rowSums(sqr(X - rowMeans(X,na.rm = TRUE)), na.rm = TRUE)/(n - 1))
}
