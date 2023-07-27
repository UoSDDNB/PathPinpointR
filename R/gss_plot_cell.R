#' GSS INDERVIDUAL CELL PLOT
#'
#' @description
#' Plots the predicted positon of a chosen cell.
#'
#' @param sample.gss A GSS_OBJECT of the sample you wish to plot
#' @param cell_idx The selected cell.
#' @param col The colour that you'd like
#' @param overlay set to TRUE if you would like this plot to overlay a previous plot.
#'
#' @return nice plot highlighting the probable position of your sample on your trajectory.
#' @export
#'
gss_cell_plot <- function(sample.gss, cell_idx = 1, col = "red", overlay = FALSE, genes_of_interest = NULL, switching_genes = NULL){

  if (!overlay) {
    plot(x = 1:100,
         y = sample.gss$cells_flat[cell_idx,],
         ylim = c(0,max(sample.gss$cells_flat[cell_idx,]) + 10),
         xlim = c(0,100),
         pch = 20,
         cex = 0.8,
         col = col,
         type = "l",
         xlab = "Pseudo-Time Index",
         ylab = "GSS Score",
         main = paste("Predicted Positions"))

    segments(which.max(sample.gss$cells_flat[cell_idx,]),
             -3.9,
             which.max(sample.gss$cells_flat[cell_idx,]),
             sample.gss$cells_flat[cell_idx,][which.max(sample.gss$cells_flat[cell_idx,])],
             lwd = 1,
             lty = 2,
             col = col)

    text(x = which.max(sample.gss$cells_flat[cell_idx,]),
         y = sample.gss$cells_flat[cell_idx,][which.max(sample.gss$cells_flat[cell_idx,])],
         labels = rownames(sample.gss$cells_flat)[cell_idx],
         col = col,
         pos = 3,
         cex = 0.69)

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
          y = sample.gss$cells_flat[cell_idx,],
          ylim = c(0,max(sample.gss$cells_flat[cell_idx,]) + 10),
          xlim = c(0,100),
          pch = 20,
          cex = 0.8,
          col = col,
          type = "l")

    segments(which.max(sample.gss$cells_flat[cell_idx,]),
             -3.9,
             which.max(sample.gss$cells_flat[cell_idx,]),
             sample.gss$cells_flat[cell_idx,][which.max(sample.gss$cells_flat[cell_idx,])],
             lwd = 1,
             lty = 2,
             col = col)

    text(x = which.max(sample.gss$cells_flat[cell_idx,]),
         y = sample.gss$cells_flat[cell_idx,][which.max(sample.gss$cells_flat[cell_idx,])],
         labels = rownames(sample.gss$cells_flat)[cell_idx],
         col = col,
         pos = 3,
         cex = 0.69)

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








