#' Performs de novo discovery of differentially reactive regions for transcriptome-wide data.
#' @param reps_A Number of replicates of group A.
#' @param reps_B Number of replicates of group B.
#' @param rl List of dataframes of reactivities for each sample.
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
#' @param processes Number of parallel processes to use.
#' @param method Character specifying either guided or de novo discovery approach.
#' @param proximity_assisted Logical, if TRUE, proximally located regions are tested together.
#' @param proximity Maximum distance between constructed regions for them to be considered proximal.
#' @param proximity_defined_length If performing a "proximity-assisted" test, minimum end-to-end length of a region to be tested.
#' @return Constructs regions, reports p-values and FDR for them.
#' @export
dStructome <- function(rl, reps_A, reps_B, batches= F, min_length = 11,
                       check_signal_strength = T, check_nucs = T, check_quality = T,
                       quality = "auto", evidence = 0, signal_strength = 0.1,
                       within_combs = NULL, between_combs= NULL, ind_regions = T, gap = 1,
                       processes = "auto", method = "denovo",
                       proximity_assisted = F, proximity = 10,
                       proximity_defined_length = 30) {

  if (is.null(names(rl)) | (length(names(rl)) != length(rl))) stop("List \'rl\' supplied to dStructome without transcript names.")

  if (processes == "auto") {
    processes = parallel::detectCores()-1
  } else {
    processes = as.numeric(processes)
    if ((processes %% 1 != 0) | (processes < 0)) stop("parameter \'processes\' supplied to dStructome must be positive integer.")
  }

  if (is.null(between_combs) | is.null(within_combs)) idcombs = getCombs(reps_A, reps_B, batches,
                                                                         between_combs, within_combs)
  if (is.null(between_combs)) between_combs = idcombs$between_combs
  if (is.null(within_combs)) within_combs = idcombs$within_combs

  if (method == "guided") {

    result <- parallel::mcmapply(function(x) {
      dStruct.guided(x, reps_A, reps_B, batches,
                     within_combs, between_combs, check_quality,
                     quality, evidence)
    }, rl, mc.cores=processes)

    res_df = data.frame(t = names(rl), pval = result)
    res_df = subset(res_df, !is.na(pval))
    res_df$FDR = p.adjust(res_df$pval, "BH")

  } else if (method == "denovo") {

    result <- parallel::mcmapply(function(x) {
      dStruct(x, reps_A, reps_B, batches, min_length,
              check_signal_strength, check_nucs, check_quality,
              quality, evidence, signal_strength,
              within_combs, between_combs, ind_regions, gap,
              get_FDR = F, proximity_assisted, proximity,
              proximity_defined_length)
    }, rl, mc.cores=processes, SIMPLIFY = F)

    if (ind_regions) {
      res_df = data.frame(t= NA, Start= NA, Stop = NA, pval= NA)
      for (i in 1:length(result)) {
        if (is.null(result[[i]])) next
        res_df = rbind(res_df, data.frame(t= names(result)[i], result[[i]]))
      }

      res_df = res_df[-1, ]
      res_df$FDR = p.adjust(res_df$pval, "BH")
    } else {
      res_df = data.frame(t= NA, Start= NA, Stop = NA)
      pvals = data.frame(t= NA, pval= NA)
      for (i in 1:length(result)) {
        if (is.null(result[[i]])) next
        pvals = rbind(pvals, data.frame(t =  names(result)[i], pval = result[[i]]$pval), stringsAsFactors= F)
        res_df = rbind(res_df, data.frame(t= names(result)[i], result[[i]]$regions))
      }

      pvals = pvals[-1, ]
      res_df = res_df[-1, ]

      pvals$FDR = p.adjust(pvals$pval, "BH")

      res_df$pval = apply(res_df, 1, function(x) subset(pvals, t == as.character(x[1]))$pval)
      res_df$FDR = apply(res_df, 1, function(x) subset(pvals, t == as.character(x[1]))$FDR)
    }
  }

  return(res_df)

}
