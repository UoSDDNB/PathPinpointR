#' Accuracy Tester
#'
#' @description
#' Use this with "sample" cells taken from the reference. It will check how close the predicted position of the cells is to the real position.
#'
#' @param reference.ppr A "PPR_OBJECT" as outputted by create_racing_lines
#' @param reference.gs The GS object containing Pseudotime.
#' @param plot Logical value. If TRUE, a histogram of inaccuracy will be plotted. Default is FALSE.
#'
#' @return If plot = FALSE, returns a data frame containing accuracy results. If plot = TRUE, returns a list containing the accuracy data frame and a plot object.
#' @export
ppr_accuracy_test <- function(reference.ppr, reference.gs, plot = TRUE) {
  # Create a data frame to store the accuracy results
  accuracy <- data.frame(
    cell_names = colnames(reference.gs), # Cell names from reference.gs
    true_position_of_cells_pseudotime = reference.gs@colData$Pseudotime, # True pseudotime values from reference.gs
    true_position_of_cells_timeIDX = NA, # Placeholder for true time indices
    predicted_position_of_cells_timeIDX = NA, # Placeholder for predicted time indices
    inaccuracy = NA # Placeholder for inaccuracy values
  )

  # Calculate the time step based on the pseudotime range
  steptime <- (max(reference.gs@colData$Pseudotime) - min(reference.gs@colData$Pseudotime)) / 100

  # Calculate the true time indices based on the pseudotime values
  accuracy$true_position_of_cells_timeIDX <- round((reference.gs@colData$Pseudotime - min(reference.gs@colData$Pseudotime)) / steptime)

  # Predict the time indices for the cells using the reference.ppr
  # Note: This assumes that the names of cells in reference.ppr$cells_flat match the cell names in reference.gs
  accuracy$predicted_position_of_cells_timeIDX[match(names(apply(reference.ppr$cells_flat, 1, which.max)), accuracy$cell_names)] <- apply(reference.ppr$cells_flat, 1, which.max)
  #may be more efficient to:
  #accuracy$predicted_position_of_cells_timeIDX <- max.col(reference.ppr$cells_flat, "first")

  # Calculate the accuracy as the absolute difference between the true and predicted time indices
  accuracy$inaccuracy <- abs(accuracy$true_position_of_cells_timeIDX - accuracy$predicted_position_of_cells_timeIDX)

  if (plot) {
    hist_plot <- hist(accuracy$inaccuracy, breaks = 100, main = "Histogram of Inaccuracy", xlab = "Inaccuracy")
    mean_inaccuracy <- mean(accuracy$inaccuracy, na.rm = TRUE)
    abline(v = mean_inaccuracy, col = "red", lwd = 1)
    text(mean_inaccuracy, max(hist_plot$counts) * 0.9, labels = paste("Mean =", round(mean_inaccuracy, 2)), adj = c(0.5, 0), col = "red")
    return(list(accuracy_data = accuracy, histogram = hist_plot))
  } else {
    return(accuracy)
  }
}
