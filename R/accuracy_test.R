#' Title Accuracy Tester
#'
#' @description
#' Use this with "sample" cells taken from the reference. It will check how close the predicted position of the cells is to the real position.
#'
#' @param gss_obj A "GSS_OBJECT" as outputted by create_racing_lines
#' @param reference.sce The SCE object which has had slingshot run on it already.
#'
#' @return A score for how accurate the results are as a data frame.
#' @export
score_gss_accuracy <- function(gss_obj, reference.sce) {

  # Create a data frame to store the accuracy results
  accuracy <- data.frame(
    cell_names = colnames(reference.sce),                                         # Cell names from reference.sce
    true_position_of_cells_pseudotime = reference.sce@colData$slingPseudotime_1,  # True pseudotime values from reference.sce
    true_position_of_cells_timeIDX = NA,                                          # Placeholder for true time indices
    predicted_position_of_cells_timeIDX = NA,                                     # Placeholder for predicted time indices
    inaccuracy = NA                                                               # Placeholder for inaccuracy values
  )

  # Calculate the time step based on the pseudotime range
  steptime <- (max(reference.sce@colData$slingPseudotime_1) - min(reference.sce@colData$slingPseudotime_1)) / 100

  # Calculate the true time indices based on the pseudotime values
  accuracy$true_position_of_cells_timeIDX <- round((reference.sce@colData$slingPseudotime_1 - min(reference.sce@colData$slingPseudotime_1)) / steptime)

  # Predict the time indices for the cells using the gss_obj
  # Note: This assumes that the names of cells in gss_obj$cells_flat match the cell names in reference.sce
  accuracy$predicted_position_of_cells_timeIDX[match(names(apply(gss_obj$cells_flat, 1, which.max)), accuracy$cell_names)] <- apply(gss_obj$cells_flat, 1, which.max)

  # Calculate the accuracy as the absolute difference between the true and predicted time indices
  accuracy$inaccuracy <- abs(accuracy$true_position_of_cells_timeIDX - accuracy$predicted_position_of_cells_timeIDX)

  return(accuracy)
}
