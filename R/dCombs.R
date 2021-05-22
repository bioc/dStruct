#' Assesses within-group or between-group variation.
#' @param rdf Data.frame of reactivities for each sample.
#' @param combs Data.frame with each column containing groupings of samples.
#' @return Nucleotide-wise d score.
#' @examples
#' #Example of a data frame with reactivities.
#' reacs <- data.frame(matrix(runif(30, 0, 10), 10, 3))
#'
#' #The columns of data frame with must indicate sample grouping and id.
#' colnames(reacs) <- c("A1", "A2", "B1")
#'
#' #Get nucleotide-wise dissimilarity scores for a set of samples.
#' dCombs(rdf = reacs, combs = data.frame(c("A1", "B1")))
#' @export
dCombs <- function(rdf, combs) {
  d <- matrix(, nrow(rdf), ncol(combs))
  for (i in 1:ncol(combs)) {
    curr_comb <- as.character(combs[, i])
    curr_dat <- rdf[, curr_comb]
    d[, i] <- calcDis(curr_dat)
  }

  d <- apply(d, 1, mean, na.rm=TRUE)
  return(d)
}
