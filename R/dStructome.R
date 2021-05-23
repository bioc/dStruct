#' @title Performs de novo or guided discovery of differentially reactive regions for transcriptome-wide data.
#'
#' @description This function provides a convenient way to call the dStruct or
#' dStructGuided functions for multiple transcripts simultaneously. By default,
#' the transcripts are processed in using multiple parallel processes if available.
#'
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
#' #Run dStruct in de novo discovery mode for all the transcripts in this data in one step.
#' dStructome(lai2019, 3, 2, batches= TRUE, min_length = 21,
#'     between_combs = data.frame(c("A3", "B1", "B2")),
#'     within_combs = data.frame(c("A1", "A2", "A3")),
#'     ind_regions = TRUE, processes = 1)
#'
#' #Load data from Wan et al., 2014
#' data(wan2014)
#'
#' #Run dStruct in guide discovery mode for all the transcript regions in this data in one step.
#' dStructome(wan2014, reps_A = 2, reps_B = 1, method = "guided", processes = 1)
#' @export
dStructome <- function(rl, reps_A, reps_B, batches= FALSE, min_length = 11,
                       check_signal_strength = TRUE, check_nucs = TRUE, check_quality = TRUE,
                       quality = "auto", evidence = 0, signal_strength = 0.1,
                       within_combs = NULL, between_combs= NULL, ind_regions = TRUE, gap = 1,
                       processes = "auto", method = "denovo",
                       proximity_assisted = FALSE, proximity = 10,
                       proximity_defined_length = 30) {

  if (is.null(names(rl)) | (length(names(rl)) != length(rl))) stop("List \'rl\' supplied to dStructome without transcript names.")

  if (processes == "auto") {
    processes <- parallel::detectCores()-1
  } else {
    processes <- as.numeric(processes)
    if ((processes %% 1 != 0) | (processes < 0)) stop("parameter \'processes\' supplied to dStructome must be positive integer.")
  }

  if (is.null(between_combs) | is.null(within_combs)) idcombs <- getCombs(reps_A, reps_B, batches,
                                                                          between_combs, within_combs)
  if (is.null(between_combs)) between_combs <- idcombs$between_combs
  if (is.null(within_combs)) within_combs <- idcombs$within_combs

  if (method == "guided") {

    result <- parallel::mcmapply(function(x) {
      dStructGuided(x, reps_A, reps_B, batches,
                    within_combs, between_combs, check_quality,
                    quality, evidence)
    }, rl, mc.cores=processes)

    pvals <- result[1, ]
    names(pvals) <- NULL
    del_d <- result[2, ]
    names(del_d) <- NULL
    res_df <- data.frame(t = names(rl), pval = pvals, del_d = del_d)
    res_df <- subset(res_df, !is.na(pval))
    res_df$FDR <- p.adjust(res_df$pval, "BH")

    rownames(res_df) <- NULL

  } else if (method == "denovo") {

    result <- parallel::mcmapply(function(x) {
      dStruct(x, reps_A, reps_B, batches, min_length,
              check_signal_strength, check_nucs, check_quality,
              quality, evidence, signal_strength,
              within_combs, between_combs, ind_regions, gap,
              get_FDR = FALSE, proximity_assisted, proximity,
              proximity_defined_length)
    }, rl, mc.cores=processes, SIMPLIFY = FALSE)

    if (ind_regions) {
      res_df <- IRanges::IRanges()
      for (i in 1:length(result)) {
        if (is.null(result[[i]])) next
        this_res <- result[[i]]
        S4Vectors::mcols(this_res) <- data.frame(S4Vectors::mcols(this_res),
                                                 t = names(result)[i])
        res_df <- c(res_df, this_res)
      }

      S4Vectors::mcols(res_df) <- data.frame(S4Vectors::mcols(res_df),
                                             FDR = p.adjust(res_df@elementMetadata$pval,
                                                            "BH"))
    } else {
      res_df <- IRanges::IRanges()
      pvals_and_del_d <- data.frame(t= NA, pval= NA, del_d = NA)
      for (i in 1:length(result)) {
        if (is.null(result[[i]])) next
        pvals_and_del_d <- rbind(pvals_and_del_d,
                                 data.frame(t =  names(result)[i], pval = result[[i]]$pval,
                                            del_d = result[[i]]$del_d), stringsAsFactors= FALSE)
        this_res <- result[[i]]$regions
        S4Vectors::mcols(this_res) <- data.frame(t = rep(names(result)[i],
                                                         length(this_res)))
        res_df <- c(res_df, this_res)
      }

      pvals_and_del_d <- pvals_and_del_d[-1, ]

      pvals_and_del_d$FDR <- p.adjust(pvals_and_del_d$pval, "BH")

      res_df_with_pvals <- IRanges::IRanges()
      for (i in unique(res_df@elementMetadata$t)) {
        this_res <- res_df[res_df@elementMetadata$t == i]
        S4Vectors::mcols(this_res) <- data.frame(t = rep(i, length(this_res)),
                                                 pval = rep(subset(pvals_and_del_d,
                                                                   t == i)$pval,
                                                            length(this_res)),
                                                 del_d = rep(subset(pvals_and_del_d,
                                                                    t == i)$del_d,
                                                             length(this_res)),
                                                 FDR = rep(subset(pvals_and_del_d,
                                                                  t == i)$FDR,
                                                           length(this_res)))
        res_df_with_pvals <- c(res_df_with_pvals, this_res)
      }

      res_df <- res_df_with_pvals
    }
  }

  return(res_df)

}
