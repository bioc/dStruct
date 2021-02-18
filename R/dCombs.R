#' Assesses within-group or between-group variation.
#' @param rdf Data.frame of reactivities for each sample.
#' @param combs Data.frame with each column containing groupings of samples.
#' @return Nucleotide-wise d score.
#' @export
dCombs <- function(rdf, combs) {
  d = matrix(, nrow(rdf), ncol(combs))
  for (i in 1:ncol(combs)) {
    curr_comb = as.character(combs[, i])
    curr_dat = rdf[, curr_comb]
    d[, i] = calcDis(curr_dat)
  }

  d = apply(d, 1, mean, na.rm=T)
  return(d)
}
