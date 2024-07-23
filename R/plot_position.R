#' PLOT POSITION
#'
#' @description
#' Plots the predicted position of your sample.
#'
#' @param sample_ppr A PPR_OBJECT of the sample you wish to plot
#' @param col The colour that you'd like
#' @param overlay TRUE if you would like this plot to overlay a previous plot.
#' @param label string that you would like to assign as the label to the line.
#' @param genes_of_interest character,
#' The names of any genes that you'd like to include the switching point of.
#' @param switching_genes a matrix containing the switching gene information,
#' as produced by GeneSwitches.
#' @param raw_values logical,
#' TRUE if you want to plot the raw values instead of percentages.
#'
#' @return a plot which highlights the probable position of your sample,
#' on your trajectory.
#' @importFrom graphics segments text lines
#' @export
#'
plot_position <- function(sample_ppr,
                          col = "red",
                          overlay = FALSE,
                          label = "sample name",
                          genes_of_interest = NULL,
                          switching_genes,
                          raw_values = FALSE) {

  sample_flat <- sample_ppr$sample_flat

  # Function to add gene segments and labels
  add_gene_segments <- function(genes, switching_genes) {
    for (gene_name in genes) {
      idx <- switching_genes[gene_name, "switch_at_timeidx"]
      segments(idx,
               -3.9,
               idx,
               -0.5,
               lwd = 1,
               lty = 2)
      text(x = idx + 3,
           y = 1,
           labels = gene_name,
           srt = -20,
           cex = 0.86,
           pos = 2)
    }
  }

  # Function to plot the sample line
  plot_sample <- function(sample_flat, col, label, overlay, raw_values) {
    x_vals <- 1:100
    y_vals <- if (raw_values) {
      sample_flat 
    } else {
      (sample_flat / max(sample_flat) * 100)
    }

    if (!overlay) {
      ylim <- if (raw_values) range(y_vals) else c(0, 110)
      plot(x = x_vals, y = y_vals, ylim = ylim, xlim = c(0, 100),
           pch = 20, cex = 0.8, col = col, type = "l",
           xlab = "Pseudo-Time Index",
           ylab = if (raw_values) "Raw ppr Score" else "ppr Score (%)",
           main = "Predicted Positions")
    } else {
      lines(x = x_vals, y = y_vals, col = col)
    }

    mid_max_idx <- which_mid_max(colSums(sample_flat))
    segments(mid_max_idx,
             -3.9,
             mid_max_idx,
             y_vals[mid_max_idx],
             lwd = 1,
             lty = 2,
             col = col)
    text(x = mid_max_idx,
         y = y_vals[mid_max_idx],
         labels = label,
         col = col,
         pos = 3)
  }

  # Main plotting
  plot_sample(sample_flat, col, label, overlay, raw_values)

  # Add gene switching points if provided
  if (!is.null(genes_of_interest) && length(genes_of_interest) > 0) {
    add_gene_segments(genes_of_interest, switching_genes)
  }
}
