#' @title Returns normalizer for reactivity vector.
#'
#' @description Assesses normalization factor for raw reactivities using the 2-8 \%
#' method. Given a reactivity profile, first, remove 2\% of the nucleotides with
#' the highest reactivities. Then, the normalization factor is the mean of reactivities
#' of the 8\% of the nucleotides with the next highest reactivities.
#'
#' @param raw.estimates A vector of raw reactivities.
#' @return The normalization factor.
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
#' normalizer(c(NA, rnorm(20, 0.5, 0.3), NA, -999))
#' @export
normalizer <- function(raw.estimates) {
  raw.estimates[which(raw.estimates == -999)] <- NA
  sorted <- raw.estimates[order(raw.estimates)]
  if (any(is.na(sorted))) {
    normalize.range <- c(round((min(which(is.na(sorted)))-1) * .9), round((min(which(is.na(sorted)))-1) * .98))
  } else {
    normalize.range <- c(round(length(sorted) * .9), round(length(sorted)* .98))
  }
  normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  return(normalizer)
}
