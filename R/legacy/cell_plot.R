#' @title Legacy Function: PPR INDIVIDUAL CELL PLOT
#'
#' @description
#' Plots the predicted position of a chosen cell.
#'
#' @param sample_ppr A ppr_OBJECT of the sample you wish to plot
#' @param cell_idx The selected cell.
#' @param col The colour that you'd like
#' @param overlay TRUE if you would like this plot to overlay a previous plot.
#' @param true_position do you want to plot the "true position of your cell"
#' @param accuracy_data dataset which contains the accuracy data
#' @param switching_genes The data which includes all of the switching genes
#' @param genes_of_interest The genes that you would like to plot
#'
#' @return plot of the probable position of your sample on your trajectory.
#' @importFrom graphics segments text lines
#' @export
#'
cell_plot <- function(sample_ppr,
                          cell_idx = 1,
                          col = "red",
                          overlay = FALSE,
                          genes_of_interest = NULL,
                          switching_genes = NULL,
                          true_position = FALSE,
                          accuracy_data = NULL) {
 usethis::deprecated("This is a legacy function and will not be maintained actively.")

  # Extract data for the selected cell
  cell_data <- sample_ppr$cells_flat[cell_idx, ]
  # Find the index of the maximum value in cell data
  max_idx <- which_mid_max(cell_data)

  # Get the true position index if accuracy data is provided
  true_pos_idx <- if (!is.null(accuracy_data)) {
    accuracy_data$true_position_of_cells_timeIDX[cell_idx]
  } else {
    NA
  }

  # Get the maximum value at the max index
  max_val <- if (!is.null(cell_data[max_idx])) {
    cell_data[max_idx]
  } else {
    NA
  }

  # Get the true name of the cell if accuracy data is provided
  true_name <- if (!is.null(accuracy_data)) {
    accuracy_data$cell_names[cell_idx]
  } else {
    NA
  }

  # If overlay is FALSE, create a new plot
  if (!overlay) {
    plot(x = 1:100,
         y = cell_data,
         ylim = c(0, max(sample_ppr$cells_flat) + 10),
         xlim = c(0, 100),
         pch = 20,
         cex = 0.8,
         col = col,
         type = "l",
         xlab = "Pseudo-Time Index",
         ylab = "PPR Score",
         main = paste("Cell Positions"))

    # Add a vertical dashed line at the max index
    segments(max_idx,
             -99999,
             max_idx,
             max_val,
             lwd = 1,
             lty = 2,
             col = col)

    # Add text to indicate the predicted position
    text(x = max_idx,
         y = max_val / 1.25,
         labels = paste(rownames(sample_ppr$cells_flat)[cell_idx],
                        "(PREDICTED)"),
         col = col,
         pos = 2,
         cex = 0.69,
         srt = 90)

    # If true position is TRUE, add the true position lines and text
    if (true_position == TRUE) {
      segments(true_pos_idx,
               -99999,
               true_pos_idx,
               99999,
               lwd = 1,
               lty = 2,
               col = col)

      text(x = true_pos_idx,
           y = max_val / 1.25,
           labels = paste(true_name, "(TRUE)"),
           col = col,
           pos = 2,
           cex = 0.69,
           srt = 90)

      segments(true_pos_idx,
               max_val,
               max_idx,
               max_val,
               lwd = 1,
               lty = 1,
               col = col)

      text(x = (max_idx - true_pos_idx) / 2 + true_pos_idx,
           y = max_val,
           labels = paste(abs(max_idx - true_pos_idx), "\nINACCURACY"),
           col = col,
           pos = 3,
           cex = 0.69)
    }

    # Plot the switching genes if any are provided
    if (length(genes_of_interest) > 0) {
      for (gene_name in genes_of_interest) {
        switch_idx <- switching_genes[gene_name, "switch_at_timeidx"]
        segments(switch_idx,
                 -3.9,
                 switch_idx,
                 -0.5,
                 lwd = 1,
                 lty = 2)

        text(x = switch_idx + 3,
             y = 1,
             labels = gene_name,
             srt = -20,
             cex = 0.86,
             pos = 2)
      }
    }

  } else {
    # If overlay is TRUE, add lines to an existing plot
    lines(x = 1:100,
          y = cell_data,
          ylim = c(0, max(cell_data) + 10),
          xlim = c(0, 100),
          pch = 20,
          cex = 0.8,
          col = col,
          type = "l")

    # Add a vertical dashed line at the max index
    segments(max_idx,
             -99999,
             max_idx,
             max_val,
             lwd = 1,
             lty = 2,
             col = col)

    # Add text to indicate the predicted position
    text(x = max_idx,
         y = max_val / 1.25,
         labels = paste(rownames(sample_ppr$cells_flat)[cell_idx],
                        "(PREDICTED)"),
         col = col,
         pos = 2,
         cex = 0.69,
         srt = 90)

    # If true position is TRUE, add the true position lines and text
    if (true_position == TRUE) {
      segments(true_pos_idx,
               -99999,
               true_pos_idx,
               99999,
               lwd = 1,
               lty = 2,
               col = col)

      text(x = true_pos_idx,
           y = max_val / 1.25,
           labels = paste(true_name, "(TRUE)"),
           col = col,
           pos = 2,
           cex = 0.69,
           srt = 90)

      segments(true_pos_idx,
               max_val,
               max_idx,
               max_val,
               lwd = 1,
               lty = 1,
               col = col)

      text(x = (max_idx - true_pos_idx) / 2 + true_pos_idx,
           y = max_val,
           labels = paste(abs(max_idx - true_pos_idx), "\nINACCURACY"),
           col = col,
           pos = 3,
           cex = 0.69)
    }

    # Plot the switching genes if any are provided
    if (length(genes_of_interest) > 0) {
      for (gene_name in genes_of_interest) {
        switch_idx <- switching_genes[gene_name, "switch_at_timeidx"]
        segments(switch_idx,
                 -3.9,
                 switch_idx,
                 -0.5,
                 lwd = 1,
                 lty = 2)

        text(x = switch_idx + 3,
             y = 1,
             labels = gene_name,
             srt = -20,
             cex = 0.86,
             pos = 2)
      }
    }
  }
}
