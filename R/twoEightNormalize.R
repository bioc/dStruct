#' @title Normalizes reactivity vector.
#'
#' @description Given a reactivity profile, first, remove 2\% of the nucleotides with
#' the highest reactivities. Then, the normalization factor is the mean of reactivities
#' of the 8\% of the nucleotides with the next highest reactivities. The raw reactivities
#' are divided by the normalization factor to get normalized reactivities. This is called
#' as 2-8 \% normalization and has been a common way to normalize data from RNA
#' structurome profiling technologies such as SHAPE-Seq, Structure-Seq, etc. (see
#' Low and Weeks, 2010, Sloma et al., 2015, and Choudhary et al., 2017).
#'
#' @param raw.estimates A vector of raw reactivities.
#' @return A vector of normalized reactivities.
#'
#' @author Krishna Choudhary
#'
#' @references
#' Low JT, Weeks KM. SHAPE-directed RNA secondary structure prediction.
#' \emph{Methods}. 2010; 52(2):150–8.
#'
#' Sloma MF, Mathews DH, Chen SJ, Burke-Aguero DH. Chapter four – improving RNA
#' secondary structure prediction with structure mapping data. In:
#' \emph{Methods in Enzymology}, vol. 553. Cambridge: Academic Press: 2015. p. 91–114.
#'
#' Choudhary K, Deng F, Aviran S. Comparative and integrative analysis of RNA
#' structural profiling data: current practices and emerging questions. \emph{Quant Biol.}
#' 2017; 5(1):3–24.
#'
#' @examples
#' twoEightNormalize(c(NA, rnorm(20, 0.5, 0.3), NA, -999))
#' @export
twoEightNormalize <- function(raw.estimates) {
  normalizer <- normalizer(raw.estimates)
  raw.estimates[which(raw.estimates != -999)] <- raw.estimates[which(raw.estimates != -999)] / normalizer
  return(raw.estimates)
}
