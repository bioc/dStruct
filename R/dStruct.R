#' Performs de novo discovery of differentially reactive regions.
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
#' @return Constructs regions, reports p-values and FDR for them.
#' @export
dStruct <- function(rdf, reps_A, reps_B, batches = F, min_length = 11,
                    check_signal_strength = T, check_nucs = T, check_quality = T,
                    quality = "auto", evidence = 0, signal_strength = 0.1,
                    within_combs = NULL, between_combs= NULL, ind_regions = T, gap = 1,
                    get_FDR = T, proximity_assisted = F, proximity = 10,
                    proximity_defined_length = 30) {

  if ((quality == "auto") & min(c(reps_A, reps_B)) != 1) quality = 0.5 else if (quality == "auto") quality = 0.2
  if (is.null(between_combs) | is.null(within_combs)) idcombs = getCombs(reps_A, reps_B, batches,
                                                                         between_combs, within_combs)
  if (is.null(between_combs)) between_combs = idcombs$between_combs
  if (is.null(within_combs)) within_combs = idcombs$within_combs

  d_within = dCombs(rdf, within_combs)
  d_between = dCombs(rdf, between_combs)

  to_test = getRegions(d_within, d_between, rdf, min_length,
                       check_signal_strength, check_nucs, check_quality,
                       quality, evidence, signal_strength)

  if (is.null(to_test) | !length(to_test)) return(NULL)

  contigs_test = getContigRegions(to_test, gap)

  if (!ind_regions & !proximity_assisted) {

    result <- tryCatch({
      wilcox.test(d_within[to_test], d_between[to_test],
                  alternative = "less", paired= T)$p.value
    }, error= function(e) {
      #Place holder for those transcripts that can't be tested due to insufficient data points.
      result = NA
    })

    result = list(regions= contigs_test, pval = result)

  } else if (!proximity_assisted) {

    pvals = c()
    for (i in 1:nrow(contigs_test)) {
      curr_res = tryCatch({
        wilcox.test(d_within[contigs_test$Start[i]:contigs_test$Stop[i]],
                    d_between[contigs_test$Start[i]:contigs_test$Stop[i]],
                    alternative = "less", paired= T)$p.value

      }, error= function(e) {
        #Place holder for those transcripts that can't be tested due to insufficient data points.
        result = NA
      })
      pvals = c(pvals, curr_res)


    }

    contigs_test = cbind(contigs_test, pval = pvals)
    if (get_FDR) contigs_test = cbind(contigs_test, FDR = p.adjust(pvals, "BH"))
    result = contigs_test
  } else {

    proximal_contigs = which(contigs_test$Start[-1]-
                               contigs_test$Stop[-nrow(contigs_test)] < proximity)+1

    if (length(proximal_contigs)) {
      proximally_tied_regs = getContigRegions(proximal_contigs)
      proximally_tied_regs$Start = proximally_tied_regs$Start-1
      proximal_contigs = unlist(apply(proximally_tied_regs, 1,
                                      function(x) x[1]:x[2]))
      non_proximal = setdiff(1:nrow(contigs_test),
                             proximal_contigs)
    } else {
      non_proximal = 1:nrow(contigs_test)
    }

    which_too_short = apply(contigs_test[non_proximal, ], 1,
                            function(x) x[2]- x[1] + 1 < proximity_defined_length)

    pvals = c()
    for (i in non_proximal) {

      curr_nucs = contigs_test$Start[i]:contigs_test$Stop[i]
      if (length(curr_nucs) < proximity_defined_length) {
        curr_res = NA
      } else {
        curr_res = tryCatch({
          wilcox.test(d_within[curr_nucs],
                      d_between[curr_nucs],
                      alternative = "less", paired= T)$p.value

        }, error= function(e) {
          #Place holder for those transcripts that can't be tested due to insufficient data points.
          result = NA
        })
      }

      pvals = c(pvals, curr_res)
    }

    if (length(proximal_contigs)) {
      proximity_defined_regs = list(Start = c(), Stop = c())
      for (i in 1:nrow(proximally_tied_regs)) {
        curr = contigs_test[proximally_tied_regs$Start[i]:proximally_tied_regs$Stop[i], ]
        curr_nucs = unlist(apply(curr, 1,
                                 function(x) x[1]:x[2]))
        curr_length = max(curr_nucs) - min(curr_nucs) + 1

        if (curr_length >= proximity_defined_length) {
          curr_res = tryCatch({
            wilcox.test(d_within[curr_nucs],
                        d_between[curr_nucs],
                        alternative = "less", paired= T)$p.value

          }, error= function(e) {
            #Place holder for those transcripts that can't be tested due to insufficient data points.
            result = NA
          })
        } else {
          curr_res = NA
        }

        # print(paste(curr$Start, collapse = ","))
        proximity_defined_regs$Start = c(proximity_defined_regs$Start,
                                         paste(curr$Start, collapse = ","))
        proximity_defined_regs$Stop = c(proximity_defined_regs$Stop,
                                        paste(curr$Stop, collapse = ","))
        pvals = c(pvals, curr_res)
      }

      contigs_test = rbind(contigs_test[non_proximal, ],
                           as.data.frame(proximity_defined_regs))
    }

    contigs_test = cbind(contigs_test, pval = pvals)
    if (get_FDR) contigs_test = cbind(contigs_test, FDR = p.adjust(pvals, "BH"))
    result = contigs_test
  }

  return(result)

}
