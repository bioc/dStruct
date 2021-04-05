#' Calculates d score.
#' @param x A numeric vector or matrix.
#' @return If input is a numeric vector, a number is returned. For a matrix, a numeric vector is returned.
#' @export
calcDis <- function(x) {
  if (class(x)== "numeric") return(2*(atan(abs(sd(x)/mean(x))))/pi )
  return(apply(x, 1, function(y) (2*(atan(abs(sd(y)/mean(y))))/pi ) ))
}
