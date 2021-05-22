#' @title Performs de novo discovery of differentially reactive regions.
#'
#' @description This function takes reactivity profiles for samples of two groups
#' as input and identifies differentially reactive regions in three steps (see
#' Choudhary et al., \emph{Genome Biology}, 2019 for details). First, it regroups
#' the samples into homogeneous and heteregenous sub-groups, which are used to
#' compute the within-group and between-group nucleotide-wise \emph{d} scores.
#' Second, smoothed between- and within-group \emph{d} score profiles are compared
#' to construct candidate differential regions. Finally, unsmoothed between- and
#' within-group \emph{d} scores are compared using the Wilcoxon signed-rank test.
#' The resulting p-values quantify the significance of difference in reactivity
#' patterns between the two input groups.
#'
#' @param reps_A Number of replicates of group A.
#' @param reps_B Number of replicates of group B.
#' @param rdf Dataframe of reactivities for each sample.
#' @param min_length Minimum length of constructed regions.
#' @param check_signal_strength Logical, if TRUE, construction of regions must be based on nucleotides that have a minimum absolute value of reactivity.
#' @param check_nucs Logical, if TRUE, constructed regions must have a minimum number of nucleotides participating in Wilcoxon signed rank test.
#' @param check_quality Logical, if TRUE, check constructed regions for quality.
#' @param batches Logical suggesting if replicates of group A and B were performed in batches and are labelled accordingly. If TRUE, a heterogeneous/homogeneous subset may not have multiple samples from the same batch.
#' @param between_combs Dataframe with each column containing groupings of replicates of groups A and B, which will be used to assess between-group variation.
#' @param within_combs Data.frame with each column containing groupings of replicates of groups A or B, which will be used to assess within-group variation.
#' @param quality Worst allowed quality for a region to be tested.
#' @param evidence Minimum evidence of increase in variation from within-group comparisons to between-group comparisons for a region to be tested.
#' @param signal_strength Threshold for minimum signal strength.
#' @param ind_regions Logical, if TRUE, test each region found in the transcript separately.
#' @param gap Integer. Join regions if they are separated by these many nucleotides.
#' @param get_FDR Logical, if FALSE, FDR is not reported.
#' @param proximity_assisted Logical, if TRUE, proximally located regions are tested together.
#' @param proximity Maximum distance between constructed regions for them to be considered proximal.
#' @param proximity_defined_length If performing a "proximity-assisted" test, minimum end-to-end length of a region to be tested.
#' @return Constructs regions, reports p-value and median difference of between-group and within-group d-scores for each region, and FDR for them.
#'
#' @author Krishna Choudhary
#'
#' @references
#' Choudhary, K., Lai, Y. H., Tran, E. J., & Aviran, S. (2019).
#' dStruct: identifying differentially reactive regions from RNA structurome
#' profiling data. \emph{Genome biology}, 20(1), 1-26.
#'
#' @examples
#' #Load data from Lai et al., 2019
#' data(lai2019)
#'
#' #Run dStruct in de novo discovery mode for a transcript with id YAL042W.
#' dStruct(rdf = lai2019[["YAL042W"]], reps_A = 3, reps_B = 2,
#'     batches = TRUE, min_length = 21,
#'     between_combs = data.frame(c("A3", "B1", "B2")),
#'     within_combs = data.frame(c("A1", "A2", "A3")),
#'     ind_regions = TRUE)
#' @export
dStruct <- function(rdf, reps_A, reps_B, batches = FALSE, min_length = 11,
                    check_signal_strength = TRUE, check_nucs = TRUE, check_quality = TRUE,
                    quality = "auto", evidence = 0, signal_strength = 0.1,
                    within_combs = NULL, between_combs= NULL, ind_regions = TRUE, gap = 1,
                    get_FDR = TRUE, proximity_assisted = FALSE, proximity = 10,
                    proximity_defined_length = 30) {

  if ((quality == "auto") & min(c(reps_A, reps_B)) != 1) quality <- 0.5 else if (quality == "auto") quality <- 0.2
  if (is.null(between_combs) | is.null(within_combs)) idcombs <- getCombs(reps_A, reps_B, batches,
                                                                         between_combs, within_combs)
  if (is.null(between_combs)) between_combs <- idcombs$between_combs
  if (is.null(within_combs)) within_combs <- idcombs$within_combs

  d_within <- dCombs(rdf, within_combs)
  d_between <- dCombs(rdf, between_combs)

  to_test <- getRegions(d_within, d_between, rdf, min_length,
                       check_signal_strength, check_nucs, check_quality,
                       quality, evidence, signal_strength)

  if (is.null(to_test) | !length(to_test)) return(NULL)

  contigs_test <- getContigRegions(to_test, gap)

  if (!ind_regions & !proximity_assisted) {

    result <- tryCatch({
      c(wilcox.test(d_within[to_test], d_between[to_test],
                    alternative = "less", paired= TRUE)$p.value,
        median(d_between[to_test] - d_within[to_test],
               na.rm = TRUE))
    }, error= function(e) {
      #Place holder for those transcripts that can't be tested due to insufficient data points.
      result <- c(NA, NA)
    })

    result <- list(regions= contigs_test, pval = result[1],
                  del_d = result[2])

  } else if (!proximity_assisted) {

    pvals <- c()
    del_d <- c()
    for (i in 1:nrow(contigs_test)) {
      curr_res <- tryCatch({
        c(wilcox.test(d_within[contigs_test$Start[i]:contigs_test$Stop[i]],
                      d_between[contigs_test$Start[i]:contigs_test$Stop[i]],
                      alternative = "less", paired= TRUE)$p.value,
          median(d_between[contigs_test$Start[i]:contigs_test$Stop[i]] -
                   d_within[contigs_test$Start[i]:contigs_test$Stop[i]],
                 na.rm = TRUE))

      }, error= function(e) {
        #Place holder for those transcripts that can't be tested due to insufficient data points.
        result <- c(NA, NA)
      })
      pvals <- c(pvals, curr_res[1])
      del_d <- c(del_d, curr_res[2])

    }

    contigs_test <- cbind(contigs_test, pval = pvals, del_d = del_d)
    if (get_FDR) contigs_test <- cbind(contigs_test, FDR = p.adjust(pvals, "BH"))
    result <- contigs_test
  } else {

    proximal_contigs <- which(contigs_test$Start[-1]-
                               contigs_test$Stop[-nrow(contigs_test)] < proximity)+1

    if (length(proximal_contigs)) {
      proximally_tied_regs <- getContigRegions(proximal_contigs)
      proximally_tied_regs$Start <- proximally_tied_regs$Start-1
      proximal_contigs <- unlist(apply(proximally_tied_regs, 1,
                                      function(x) x[1]:x[2]))
      non_proximal <- setdiff(1:nrow(contigs_test),
                             proximal_contigs)
    } else {
      non_proximal <- 1:nrow(contigs_test)
    }

    which_too_short <- apply(contigs_test[non_proximal, ], 1,
                            function(x) x[2]- x[1] + 1 < proximity_defined_length)

    pvals <- c()
    del_d <- c()
    for (i in non_proximal) {

      curr_nucs <- contigs_test$Start[i]:contigs_test$Stop[i]
      if (length(curr_nucs) < proximity_defined_length) {
        curr_res <- c(NA, NA)
      } else {
        curr_res <- tryCatch({
          c(wilcox.test(d_within[curr_nucs],
                        d_between[curr_nucs],
                        alternative = "less", paired= TRUE)$p.value,
            median(d_between[curr_nucs] - d_within[curr_nucs],
                   na.rm = TRUE))

        }, error= function(e) {
          #Place holder for those transcripts that can't be tested due to insufficient data points.
          result <- c(NA, NA)
        })
      }

      pvals <- c(pvals, curr_res[1])
      del_d <- c(del_d, curr_res[2])
    }

    if (length(proximal_contigs)) {
      proximity_defined_regs <- list(Start = c(), Stop = c())
      for (i in 1:nrow(proximally_tied_regs)) {
        curr <- contigs_test[proximally_tied_regs$Start[i]:proximally_tied_regs$Stop[i], ]
        curr_nucs <- unlist(apply(curr, 1,
                                 function(x) x[1]:x[2]))
        curr_length <- max(curr_nucs) - min(curr_nucs) + 1

        if (curr_length >= proximity_defined_length) {
          curr_res <- tryCatch({
            c(wilcox.test(d_within[curr_nucs],
                        d_between[curr_nucs],
                        alternative = "less", paired= TRUE)$p.value,
            median(d_between[curr_nucs] - d_within[curr_nucs],
                   na.rm = TRUE))

          }, error= function(e) {
            #Place holder for those transcripts that can't be tested due to insufficient data points.
            result <- c(NA, NA)
          })
        } else {
          curr_res <- c(NA, NA)
        }

        proximity_defined_regs$Start <- c(proximity_defined_regs$Start,
                                         paste(curr$Start, collapse = ","))
        proximity_defined_regs$Stop <- c(proximity_defined_regs$Stop,
                                        paste(curr$Stop, collapse = ","))
        pvals <- c(pvals, curr_res[1])
        del_d <- c(del_d, curr_res[2])
      }

      contigs_test <- rbind(contigs_test[non_proximal, ],
                           as.data.frame(proximity_defined_regs))
    }

    contigs_test <- cbind(contigs_test, pval = pvals, del_d = del_d)
    if (get_FDR) contigs_test <- cbind(contigs_test, FDR = p.adjust(pvals, "BH"))
    result <- contigs_test
  }

  return(result)

}
