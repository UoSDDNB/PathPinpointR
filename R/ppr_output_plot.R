#' PPR OUTPUT PLOT
#'
#' @description
#' Plots the predicted positon of your sample.
#'
#' @param sample_ppr A PPR_OBJECT of the sample you wish to plot
#' @param col The colour that you'd like
#' @param overlay set to TRUE if you would like this plot to overlay a previous plot.
#' @param label string that you would like to assign as the label to the line.
#' @param genes_of_interest The names of any genes that you'd like to include the switching point of.
#' @param switching_genes a matrix containing the switching gene information as produced by GeneSwitches.
#'
#' @return nice plot highlighting the probable position of your sample on your trajectory.
#' @importFrom graphics segments text lines
#' @export
#'
ppr_output_plot <- function(sample_ppr, col = "red", overlay = FALSE, label = "sample name", genes_of_interest = NULL, switching_genes){

 if (!overlay) {
  plot(x = 1:100,
       y = (sample_ppr$sample_flat/max(sample_ppr$sample_flat)*100),
       ylim = c(0,110),
       xlim = c(0,100),
       pch = 20,
       cex = 0.8,
       col = col,
       type = "l",
       xlab = "Pseudo-Time Index",
       ylab = "ppr Score",
       main = paste("Predicted Positions"))

  segments(which_mid_max(colSums(sample_ppr$sample_flat)),
           -3.9,
           which_mid_max(colSums(sample_ppr$sample_flat)),
           (sample_ppr$sample_flat[which_mid_max(colSums(sample_ppr$sample_flat))]/max(sample_ppr$sample_flat)*100),
           lwd = 1,
           lty = 2,
           col = col)

  text(x = which_mid_max(colSums(sample_ppr$sample_flat)),
       y = (sample_ppr$sample_flat[which_mid_max(colSums(sample_ppr$sample_flat))]/max(sample_ppr$sample_flat)*100),
       labels = label,
       col = col,
       pos = 3)

  if (length(genes_of_interest) > 0) {
    for (gene_name in genes_of_interest) {
      segments(switching_genes[gene_name,"switch_at_timeidx"],
               -3.9,
               switching_genes[gene_name,"switch_at_timeidx"],
               -0.5,
               lwd = 1,
               lty = 2)

      text(x = switching_genes[gene_name,"switch_at_timeidx"] + 3,
           y = 1,
           labels = gene_name,
           srt = -20,
           cex = 0.86,
           pos = 2)

    }
  }



  } else {
    lines(x = 1:100,
         y = (sample_ppr$sample_flat/max(sample_ppr$sample_flat)*100),
         ylim = c(0,110),
         xlim = c(0,100),
         pch = 20,
         cex = 0.8,
         col = col,
         type = "l")

    segments(which_mid_max(colSums(sample_ppr$sample_flat)),
             -3.9,
             which_mid_max(colSums(sample_ppr$sample_flat)),
             (sample_ppr$sample_flat[which_mid_max(colSums(sample_ppr$sample_flat))]/max(sample_ppr$sample_flat)*100),
             lwd = 1,
             lty = 2,
             col = col)

    text(x = which_mid_max(colSums(sample_ppr$sample_flat)),
         y = (sample_ppr$sample_flat[which_mid_max(colSums(sample_ppr$sample_flat))]/max(sample_ppr$sample_flat)*100),
         labels = label,
         col = col,
         pos = 3)

    if (length(genes_of_interest) > 0) {
      for (gene_name in genes_of_interest) {
        segments(switching_genes[gene_name,"switch_at_timeidx"],
                 -3.9,
                 switching_genes[gene_name,"switch_at_timeidx"],
                 -0.5,
                 lwd = 1,
                 lty = 2)

        text(x = switching_genes[gene_name,"switch_at_timeidx"] + 3,
             y = 1,
             labels = gene_name,
             srt = -20,
             cex = 0.86,
             pos = 2)

   }
  }
 }
}








