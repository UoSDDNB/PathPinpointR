#' Accuracy Tester
#'
#' @description
#' Use this with "sample" cells taken from the reference. It will check how close the predicted position of the cells is to the real position.
#'
#' @param reference.pppr A "GSS_OBJECT" as outputted by create_racing_lines
#' @param reference.gs The GS object containing Pseudotime.
#'
#' @return A score for how accurate the results are as a data frame.
#' @export
pppr_accuracy_test <- function(reference.pppr, reference.gs) {

  # Create a data frame to store the accuracy results
  accuracy <- data.frame(
    cell_names = colnames(reference.gs),                                          # Cell names from reference.gs
    true_position_of_cells_pseudotime = reference.gs@colData$Pseudotime,          # True pseudotime values from reference.gs
    true_position_of_cells_timeIDX = NA,                                          # Placeholder for true time indices
    predicted_position_of_cells_timeIDX = NA,                                     # Placeholder for predicted time indices
    inaccuracy = NA                                                               # Placeholder for inaccuracy values
  )

  # Calculate the time step based on the pseudotime range
  steptime <- (max(reference.gs@colData$Pseudotime) - min(reference.gs@colData$Pseudotime)) / 100

  # Calculate the true time indices based on the pseudotime values
  accuracy$true_position_of_cells_timeIDX <- round((reference.gs@colData$Pseudotime - min(reference.gs@colData$Pseudotime)) / steptime)

  # Predict the time indices for the cells using the reference.pppr
  # Note: This assumes that the names of cells in reference.pppr$cells_flat match the cell names in reference.gs
  accuracy$predicted_position_of_cells_timeIDX[match(names(apply(reference.pppr$cells_flat, 1, which.max)), accuracy$cell_names)] <- apply(reference.pppr$cells_flat, 1, which.max)

  # Calculate the accuracy as the absolute difference between the true and predicted time indices
  accuracy$inaccuracy <- abs(accuracy$true_position_of_cells_timeIDX - accuracy$predicted_position_of_cells_timeIDX)

  return(accuracy)
}
