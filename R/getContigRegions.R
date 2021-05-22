#' @title Identifies contiguous regions from a list of nucleotide indices.
#'
#' @description Given a sequence of nucleotide indices, this function returns integer ranges covered by the indices. There is an option to merge ranges if they are separated by less than a user-specified distance.
#'
#' @param x A vector of integers.
#' @param gap Allowed gap to merge regions.
#' @return Dataframe storing start and stop sites of continguous regions.
#'
#' @author Krishna Choudhary
#'
#' @references
#' Choudhary, K., Lai, Y. H., Tran, E. J., & Aviran, S. (2019).
#' dStruct: identifying differentially reactive regions from RNA
#' structurome profiling data. \emph{Genome biology}, 20(1), 1-26.
#'
#' @examples
#' #Convert an integer vector of nucleotide positions to a data frame containing the coordinates of contiguous regions.
#' nucleotide_positions <- c(1, 3, 2, 8, 4:7, 11:20)
#' getContigRegions(nucleotide_positions)
#'
#' #Merge regions if their end points are within 3 nt of each other.
#' getContigRegions(nucleotide_positions, gap = 3)
#' @export
getContigRegions <- function(x, gap = 0) {
  if ((gap %% 1 != 0) | (gap < 0)) stop("parameter \'gap\' supplied to getContigRegions must be positive integer.")
  x <- sort(x)
  return(data.frame(Start= c(x[1], x[which(diff(x) > 1 + gap) +1]),
                    Stop = c(x[which(diff(x) > 1 + gap)], tail(x, 1))))
}
