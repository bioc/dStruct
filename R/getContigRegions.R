#' @title Identifies contiguous regions from a list of nucleotide indices.
#'
#' @description Given a sequence of nucleotide indices, this function returns integer ranges covered by the indices. There is an option to merge ranges if they are separated by less than a user-specified distance.
#'
#' @param x A vector of integers.
#' @param gap Include gaps in the ranges if they are shorter than or equal to this length.
#' @return IRanges object storing start and end sites of continguous regions.
#'
#' @author Krishna Choudhary
#'
#' @examples
#' #Convert an integer vector of nucleotide positions to an IRanges object containing the coordinates of contiguous regions.
#' nucleotide_positions <- c(1, 3, 2, 8, 4:7, 11:20)
#' getContigRegions(nucleotide_positions)
#'
#' #Merge regions if their end points are within 3 nt of each other.
#' getContigRegions(nucleotide_positions, gap = 3)
#' @export
getContigRegions <- function(x, gap = 0) {
  if ((gap %% 1 != 0) | (gap < 0)) stop("parameter \'gap\' supplied to getContigRegions must be positive integer.")
  x <- sort(x)
  x <- replace(logical(max(x)), x, TRUE)
  ir <- IRanges::IRanges(x)

  if (gap) {
    gir <- IRanges::gaps(ir)
    gir <- gir[IRanges::width(gir) <= gap]
    ir <- IRanges::reduce(c(ir, gir))
  }

  return(ir)
}
