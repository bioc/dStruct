#' Identifies subgroupings of replicates for assessing within-group and between-group variation.
#' @param reps_A Number of replicates of group A.
#' @param reps_B Number of replicates of group B.
#' @param batches Logical suggesting if replicates of group A and B were performed in batches and are labelled accordingly. If TRUE, a heterogeneous/homogeneous subset may not have multiple samples from the same batch.
#' @param between_combs Dataframe with each column containing groupings of replicates of groups A and B, which will be used to assess between-group variation.
#' @param within_combs Data.frame with each column containing groupings of replicates of groups A or B, which will be used to assess within-group variation.
#' @return List of two dataframes, containing groupings for within-group and between-group variation.
#' @export
getCombs <- function(reps_A, reps_B, batches = FALSE, between_combs= NULL, within_combs= NULL) {
  if ((is.null(within_combs) & !is.null(between_combs)) | (is.null(between_combs) & !is.null(within_combs))) stop("Homogeneous and heterogeneous sets should either both be null or both be specified with equal number of samples in each set.")
  if (!is.null(within_combs)) if (nrow(between_combs)!=nrow(within_combs)) stop("Heterogeneous and homogeneous sets should have equal number of samples.")

  set_membership = if (is.null(within_combs)) max(c(reps_A, reps_B)) else nrow(within_combs)
  while ( (set_membership/2 > min(c(reps_A, reps_B))) | (set_membership %% 2 != 0) ) {

    if (set_membership  == 3) break

    set_membership = set_membership - 1

  }

  all_combs = combn(c(paste0("A", 1:reps_A), paste0("B", 1:reps_B)),
                    set_membership)


  if (batches) {
    invalid_due_to_batch = apply(all_combs, 2, function(x) {
      curr_reps = mapply(function(y) y[2], strsplit(x, split= ""))
      return(length(unique(curr_reps)) < length(curr_reps))
    })

    all_combs = all_combs[, !invalid_due_to_batch] #Removes sets that have multiple samples from the same batch.
  }

  if (is.null(within_combs)) {
    if (min(c(reps_A, reps_B)) >= set_membership) {
      within_combs = cbind(combn(paste0("A", 1:reps_A), set_membership),
                           combn(paste0("B", 1:reps_B), set_membership))
    } else if (reps_A >= set_membership) {
      within_combs = cbind(combn(paste0("A", 1:reps_A), set_membership))
    } else if (reps_B >= set_membership) {
      within_combs = cbind(combn(paste0("B", 1:reps_B), set_membership))
    }

  }


  if (is.null(between_combs)) {
    between_combs_not_in = which(apply(all_combs, 2, paste0, collapse= "") %in% apply(within_combs, 2, paste0, collapse= ""))
    between_combs = data.frame(all_combs[, -between_combs_not_in])
    deg_of_heterogeneity = apply(between_combs, 2, function(x) {
      curr_reps = mapply(function(y) y[1], strsplit(x, split= ""))
      gA = length(which(curr_reps == "A"))
      gB = length(which(curr_reps == "B"))
      return(gA*gB/(gA + gB))
    })

    to_keep = which(deg_of_heterogeneity == max(deg_of_heterogeneity))
    between_combs = data.frame(between_combs[, to_keep])
  }

  return(list(between_combs= between_combs, within_combs= within_combs))
}
