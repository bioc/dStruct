#' Returns normalizer for reactivity vector.
#' @description Assesses normalization factor for raw reactivities using 2-8 \% method.
#' @param raw.estimates A vector of raw reactivities.
#' @return The normalization factor.
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
