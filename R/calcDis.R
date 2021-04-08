#' Calculates d score.
#' @param x A numeric vector or matrix.
#' @return If input is a numeric vector, a number is returned. For a matrix, a numeric vector is returned.
#' @examples
#' #Lower standard deviation of reactivites results in lower d-score.
#' calcDis(rnorm(10, 1, 0.2))
#' calcDis(rnorm(10, 1, 0.6))
#' @export
calcDis <- function(x) {
  dScore <- function(y) 2*(atan(abs(sd(y)/mean(y))))/pi
  if (is.numeric(x)) return(dScore(x))
  return(apply(x, 1, function(y) dScore(y)))
}
