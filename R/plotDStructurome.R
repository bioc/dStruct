#' @title Plots differentially reactive regions.
#'
#' @description Given the table of results from dStruct or dStructGuided and the
#' corresponding lists with reactivity scores for all transcripts, this function saves
#' a PDF file with detailed visualizations of reactivities for all differential regions.
#'
#' @param rl List of dataframes of reactivities for each sample.
#' @param diff_regions Dataframe of regions with significance of differentially reactivity.
#' @param outfile The name for pdf file which will be saved.
#' @param fdr FDR threshold for plotted regions.
#' @param ylim Y-axis limits for plots.
#' @param del_d_cutoff Minimum effect size for plotted regions specified in terms of median difference of the between-group and within-group d-scores.
#' @return Saves a PDF for all differentially reactive regions. Returns NULL.
#'
#' @author Krishna Choudhary
#'
#' @references
#' Choudhary, K., Lai, Y. H., Tran, E. J., & Aviran, S. (2019).
#' dStruct: identifying differentially reactive regions from RNA
#' structurome profiling data. \emph{Genome biology}, 20(1), 1-26.
#'
#' @examples
#' #Load data from Lai et al., 2019
#' data(lai2019)
#'
#' #Run dStruct in de novo discovery mode for all the transcripts in this data in one step.
#' res <- dStructome(lai2019, 3, 2, batches= TRUE, min_length = 21,
#'     between_combs = data.frame(c("A3", "B1", "B2")),
#'     within_combs = data.frame(c("A1", "A2", "A3")),
#'     ind_regions = TRUE, processes = 1)
#'
#' #Plot the significant results and save to a PDF file.
#' plotDStructurome(rl = lai2019,
#'     diff_regions = res,
#'     outfile = "significantly_differential_regions",
#'     fdr = 0.05,
#'     ylim = c(-0.05, 3))
#' @export
plotDStructurome <- function(rl, diff_regions, outfile, fdr = 0.05, ylim = c(-0.05, 3),
                              del_d_cutoff = 0.01) {
  diff_regions <- data.frame(diff_regions)
  diff_regions <- subset(diff_regions, FDR < fdr, del_d > del_d_cutoff)
  diff_t <- unique(diff_regions$t)

  pdf(paste0(outfile, ".pdf"), width=7,height=6)
  for (i in 1:length(diff_t)) {
    curr_t <- as.character(diff_t[i])
    curr_regs <- subset(diff_regions, t == curr_t)
    curr_df <- rl[[curr_t]]

    print(ggplot2::ggplot(data=data.frame(x=0,y=0), ggplot2::aes(x=x, y=y)) +
            ggplot2::annotate("text", x = 4, y = 25,
                              label = curr_t) +
            ggplot2::theme_classic()+
            ggplot2::theme(axis.line=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks=ggplot2::element_blank(),
                           axis.title.x=ggplot2::element_blank(),
                           axis.title.y=ggplot2::element_blank(),
                           legend.position="none",
                           panel.background=ggplot2::element_blank(),
                           panel.border=ggplot2::element_blank(),
                           panel.grid.major=ggplot2::element_blank(),
                           panel.grid.minor=ggplot2::element_blank(),
                           plot.background=ggplot2::element_blank())
    )

    dat <- data.frame(curr_df, n = 1:nrow(curr_df))
    dat <- reshape2::melt(dat, id.vars = "n")

    print(ggplot2::ggplot(dat, ggplot2::aes(x= n, y= value)) +
            ggplot2::geom_bar(stat="identity") +ggplot2::facet_grid(variable~.)+
            ggplot2::coord_cartesian(ylim=ylim) +
            ggplot2::geom_rect(ggplot2::aes(NULL, NULL, xmin=start-0.5, xmax=end+0.5),
                               ymin= -Inf, ymax= Inf, data= curr_regs, fill= "red",
                               color= NA, alpha= 0.3) +
            ggplot2::theme(strip.text = ggplot2::element_text(size = 6),
                           legend.position = "none") +
            ggplot2::ggtitle(curr_t)+
            ggplot2::xlab("Nucleotide") + ggplot2::ylab("Reactivity"))


    for (r in 1:nrow(curr_regs)) {
      dat <- data.frame(curr_df[curr_regs$start[r]:curr_regs$end[r], ],
                       n = curr_regs$start[r]:curr_regs$end[r])
      dat <- reshape2::melt(dat, id.vars = "n")

      print(ggplot2::ggplot(dat, ggplot2::aes(x= n, y= value)) +
              ggplot2::geom_bar(stat="identity") +ggplot2::facet_grid(variable~.)+
              ggplot2::coord_cartesian(ylim=ylim) +
              ggplot2::theme(strip.text = ggplot2::element_text(size = 6),
                             legend.position = "none") +
              ggplot2::ggtitle(curr_t)+
              ggplot2::xlab("Nucleotide") + ggplot2::ylab("Reactivity"))
    }

  }
  dev.off()
  return()
}

