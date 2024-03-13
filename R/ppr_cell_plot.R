#' PPR INDERVIDUAL CELL PLOT
#'
#' @description
#' Plots the predicted positon of a chosen cell.
#'
#' @param sample.ppr A ppr_OBJECT of the sample you wish to plot
#' @param cell_idx The selected cell.
#' @param col The colour that you'd like
#' @param overlay set to TRUE if you would like this plot to overlay a previous plot.
#' @param true_position do you wan to plot the "true position of your cell"
#' @param accuracy_data dataset which contains the accuracy data
#' @param switching_genes The data which includes all of the switching genes
#' @param genes_of_interest The genes that you would like to plot
#'
#' @return nice plot highlighting the probable position of your sample on your trajectory.
#' @importFrom graphics segments text lines
#' @export
#'
ppr_cell_plot <- function(sample.ppr, cell_idx = 1, col = "red", overlay = FALSE, genes_of_interest = NULL, switching_genes = NULL, true_position = FALSE, accuracy_data = NULL) {
  if (!overlay) {
    plot(x = 1:100,
         y = sample.ppr$cells_flat[cell_idx,],
         ylim = c(0,max(sample.ppr$cells_flat) + 10),
         xlim = c(0,100),
         pch = 20,
         cex = 0.8,
         col = col,
         type = "l",
         xlab = "Pseudo-Time Index",
         ylab = "GSS Score",
         main = paste("Cell Positions"))

    segments(which.max(sample.ppr$cells_flat[cell_idx,]),
             -99999,
             which.max(sample.ppr$cells_flat[cell_idx,]),
             sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])],
             lwd = 1,
             lty = 2,
             col = col)

    text(x = which.max(sample.ppr$cells_flat[cell_idx,]),
         y = (sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])])/1.25,
         labels = paste(rownames(sample.ppr$cells_flat)[cell_idx],"(PREDICTED)"),
         col = col,
         pos = 2,
         cex = 0.69,
         srt = 90)

    if(true_position == TRUE) {
      segments(accuracy_data$true_position_of_cells_timeIDX[cell_idx],
               -99999,
               accuracy_data$true_position_of_cells_timeIDX[cell_idx],
               99999,
               lwd = 1,
               lty = 2,
               col = col)

      text(x = accuracy_data$true_position_of_cells_timeIDX[cell_idx],
           y = (sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])])/1.25,
           labels = paste(accuracy_data$cell_names[cell_idx],"(TRUE)"),
           col = col,
           pos = 2,    # Set the position to 2 (Left)
           cex = 0.69,
           srt = 90)   # Set the string rotation to 90 degrees

      segments(accuracy_data$true_position_of_cells_timeIDX[cell_idx],
               sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])],
               which.max(sample.ppr$cells_flat[cell_idx,]),
               sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])],
               lwd = 1,
               lty = 1,
               col = col)

      text(x = (which.max(sample.ppr$cells_flat[cell_idx,]) - accuracy_data$true_position_of_cells_timeIDX[cell_idx] )/2 +accuracy_data$true_position_of_cells_timeIDX[cell_idx],
           y = sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])],
           labels = paste(abs(which.max(sample.ppr$cells_flat[cell_idx,]) - accuracy_data$true_position_of_cells_timeIDX[cell_idx]), "\nINACCURACY"),
           col = col,
           pos = 3,    # Set the position to 3 (TOP)
           cex = 0.69)

    }


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
          y = sample.ppr$cells_flat[cell_idx,],
          ylim = c(0,max(sample.ppr$cells_flat[cell_idx,]) + 10),
          xlim = c(0,100),
          pch = 20,
          cex = 0.8,
          col = col,
          type = "l")

    segments(which.max(sample.ppr$cells_flat[cell_idx,]),
             -99999,
             which.max(sample.ppr$cells_flat[cell_idx,]),
             sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])],
             lwd = 1,
             lty = 2,
             col = col)

    text(x = which.max(sample.ppr$cells_flat[cell_idx,]),
         y = (sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])])/1.25,
         labels = paste(rownames(sample.ppr$cells_flat)[cell_idx],"(PREDICTED)"),
         col = col,
         pos = 2,
         cex = 0.69,
         srt = 90)

    if(true_position == TRUE) {
      segments(accuracy_data$true_position_of_cells_timeIDX[cell_idx],
               -99999,
               accuracy_data$true_position_of_cells_timeIDX[cell_idx],
               99999,
               lwd = 1,
               lty = 2,
               col = col)

      text(x = accuracy_data$true_position_of_cells_timeIDX[cell_idx],
           y = (sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])])/1.25,
           labels = paste(accuracy_data$cell_names[cell_idx],"(TRUE)"),
           col = col,
           pos = 2,    # Set the position to 2 (Left)
           cex = 0.69,
           srt = 90)   # Set the string rotation to 90 degrees

      segments(accuracy_data$true_position_of_cells_timeIDX[cell_idx],
               sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])],
               which.max(sample.ppr$cells_flat[cell_idx,]),
               sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])],
               lwd = 1,
               lty = 1,
               col = col)

      text(x = (which.max(sample.ppr$cells_flat[cell_idx,]) - accuracy_data$true_position_of_cells_timeIDX[cell_idx] )/2 +accuracy_data$true_position_of_cells_timeIDX[cell_idx],
           y = sample.ppr$cells_flat[cell_idx,][which.max(sample.ppr$cells_flat[cell_idx,])],
           labels = paste(abs(which.max(sample.ppr$cells_flat[cell_idx,]) - accuracy_data$true_position_of_cells_timeIDX[cell_idx]), "\nINACCURACY"),
           col = col,
           pos = 3,    # Set the position to 3 (TOP)
           cex = 0.69)

    }

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








