#' Normalizes reactivity vector.
#' @description Normalizes raw reactivities using 2-8 \% method.
#' @param raw.estimates A vector of raw reactivities.
#' @return A vector of normalized reactivities.
#' @examples
#' twoEightNormalize(c(NA, rnorm(20, 0.5, 0.3), NA, -999))
#' @export
twoEightNormalize <- function(raw.estimates) {
  normalizer <- normalizer(raw.estimates)
  raw.estimates[which(raw.estimates != -999)] <- raw.estimates[which(raw.estimates != -999)] / normalizer
  return(raw.estimates)
}
