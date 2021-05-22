#' @title Constructs potential differentially reactive regions.
#'
#' @description This function takes between- and within-group \emph{d} scores for a
#' transcript as input and identifies regions where the former is generally larger.
#' Regions that pass minimum quality and minimum signal criteria are returned.
#'
#' @param d_within Nucleotide-wise d score for within-group variation.
#' @param d_spec Nucleotide-wise d score for between-group variation.
#' @param rdf Dataframe of reactivities for each sample.
#' @param min_length Minimum length of constructed regions.
#' @param check_signal_strength Logical, if TRUE, construction of regions must be based on nucleotides that have a minimum absolute value of reactivity.
#' @param check_nucs Logical, if TRUE, constructed regions must have a minimum number of nucleotides participating in Wilcoxon signed rank test.
#' @param check_quality Logical, if TRUE, check constructed regions for quality.
#' @param quality Worst allowed quality for a region to be tested.
#' @param evidence Minimum evidence of increase in variation from within-group comparisons to between-group comparisons for a region to be tested.
#' @param signal_strength Threshold for minimum signal strength.
#' @return Integer vector of nucleotides that constitute potential differentially reactive regions.
#'
#' @author Krishna Choudhary
#'
#' @references
#' Choudhary, K., Lai, Y. H., Tran, E. J., & Aviran, S. (2019).
#' dStruct: identifying differentially reactive regions from RNA
#' structurome profiling data. \emph{Genome biology}, 20(1), 1-26.
#'
#' @export
getRegions <- function(d_within, d_spec, rdf, min_length= 11,
                       check_signal_strength = TRUE, check_nucs = TRUE, check_quality = TRUE,
                       quality = 0.5, evidence = 0, signal_strength = 0.1) {

  if (check_signal_strength) {
    insufficient_signal <- apply(rdf, 1, function(x) all(abs(x) < signal_strength, na.rm=TRUE))
    d_within[insufficient_signal] <- NA
    d_spec[insufficient_signal] <- NA
  }

  #Set min_length to next odd integer.
  min_length <- if (min_length %% 2 == 0) min_length +1 else min_length

  to_test <- c()

  #How many nucleotides with information must a constructed region have.
  min_nucs <- min(min_length, 6)-1

  smooth_evidence <- zoo::rollapply(d_spec - d_within, width = min_length,
                                   FUN= function(x) mean(x, na.rm=TRUE))

  where_evidence <- smooth_evidence > evidence
  smooth_evidence[where_evidence] <- 1
  smooth_evidence[!where_evidence] <- 0

  evidence_rle <- rle(smooth_evidence)
  consider_regs <- which((evidence_rle$values == 1) & (evidence_rle$lengths >= min_length))

  for (i in consider_regs) {
    check_n <- check_nucs
    check_q <- check_quality

    if (i == 1) {
      start <- which(!is.na(d_within - d_spec))[1] + ((min_length-1)/2)
      end <- evidence_rle$lengths[1]+((min_length-1)/2)
    } else {
      start <- sum(evidence_rle$lengths[1:(i-1)]) + ((min_length-1)/2) + 1
      end <- sum(evidence_rle$lengths[1:(i)]) + ((min_length-1)/2)
    }

    #Trim nucleotides from end with very low reactivities, only if regions detected are not too short.
    if (end - start + 1 >= 11) {
      trim_reac <- apply(round(rdf[start:end, ], 2), 1,
                        function(x) all(abs(x) <= signal_strength, na.rm= TRUE))
      if (all(trim_reac)) next

      start <- start + which(!trim_reac)[1] - 1
      end <- end - length(trim_reac) + tail(which(!trim_reac), 1)

      trim_reac <- is.na(d_within[start:end]-d_spec[start:end])
      if (all(trim_reac)) next
      start <- start + which(!trim_reac)[1] - 1
      end <- end - length(trim_reac) + tail(which(!trim_reac), 1)
    }

    if (end - start + 1 < min_length) next


    if (check_n) {
      check_n <- sum(!is.na(c(d_within[start:end] - d_spec[start:end]))) > min_nucs

      #check_q can be NA if for example, reactivity were all 0s in a region.
      if (is.na(check_n)) check_n <- FALSE
    } else check_n <- TRUE

    if (check_n) curr_test <- start:end

    if (check_q) {
      check_q <- mean(d_within[start:end], na.rm= TRUE) <= quality

      #check_q can be NA if for example, reactivity were all 0s in a region.
      if (!is.na(check_q)) {
        #If the quality is not good for entire region, search for smaller regions that are good quality.
        if (!check_q) {
          wi_d <- zoo::rollapply(d_within[start:end], width = min_length,
                                FUN= function(x) mean(x), align= "left")
          if (any(wi_d <= quality, na.rm= TRUE))  {
            which_good <- start + which(wi_d <= quality) -1
            curr_test <- unique(unlist(purrr::map(which_good, function(x) x:(x+min_length-1))))
            check_q <- TRUE
          }
        }
      } else check_q <- FALSE
    } else check_q <- TRUE

    if (check_n & check_q) to_test <- c(to_test, start:end)

  }

  return(to_test)
}
