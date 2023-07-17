#' GSS OUTPUT PLOT
#'
#' @description
#' Plots the predicted positon of your sample.
#'
#' @param sample.gss A GSS_OBJECT of the sample you wish to plot
#' @param col The colour that you'd like
#' @param overlay set to TRUE if you would like this plot to overlay a previous plot.
#'
#' @return nice plot highlighting the probable position of your sample on your trajectory.
#' @export
#'
gss_output_plot <- function(sample.gss, col = "red", overlay = FALSE, label = "sample name", genes_of_interest = c("ENSG00000248605","ENSG00000272168")){

 if (!overlay) {
  plot(x = 1:100,
       y = (sample.gss$sample_flat/max(sample.gss$sample_flat)*100),
       ylim = c(0,110),
       xlim = c(0,100),
       pch = 20,
       cex = 0.8,
       col = col,
       type = "p",
       xlab = "Pseudo-Time Index",
       ylab = "GSS Score",
       main = paste("Predicted Positions"))

  segments(which.max(colSums(sample.gss$sample_flat)),
           -3.9,
           which.max(colSums(sample.gss$sample_flat)),
           (sample.gss$sample_flat[which.max(colSums(sample.gss$sample_flat))]/max(sample.gss$sample_flat)*100),
           lwd = 1,
           lty = 2,
           col = col)

  text(x = which.max(colSums(sample.gss$sample_flat)),
       y = (sample.gss$sample_flat[which.max(colSums(sample.gss$sample_flat))]/max(sample.gss$sample_flat)*100),
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
         y = (sample.gss$sample_flat/max(sample.gss$sample_flat)*100),
         ylim = c(0,110),
         xlim = c(0,100),
         pch = 20,
         cex = 0.8,
         col = col,
         type = "p")

    segments(which.max(colSums(sample.gss$sample_flat)),
             -3.9,
             which.max(colSums(sample.gss$sample_flat)),
             (sample.gss$sample_flat[which.max(colSums(sample.gss$sample_flat))]/max(sample.gss$sample_flat)*100),
             lwd = 1,
             lty = 2,
             col = col)

    text(x = which.max(colSums(sample.gss$sample_flat)),
         y = (sample.gss$sample_flat[which.max(colSums(sample.gss$sample_flat))]/max(sample.gss$sample_flat)*100),
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








