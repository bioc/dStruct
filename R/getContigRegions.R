#' Identifies contiguous regions from a list of nucleotide indices.
#' @param x A vector of integers.
#' @param gap Allowed gap to merge regions.
#' @return Dataframe storing start and stop sites of continguous regions.
#' @export
getContigRegions <- function(x, gap = 0) {
  if ((gap %% 1 != 0) | (gap < 0)) stop("parameter \'gap\' supplied to getContigRegions must be positive integer.")
  x = sort(x)
  return(data.frame(Start= c(x[1], x[which(diff(x) > 1 + gap) +1]),
                    Stop = c(x[which(diff(x) > 1 + gap)], tail(x, 1))))
}
