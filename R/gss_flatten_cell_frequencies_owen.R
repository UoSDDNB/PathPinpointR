#' Title flatten_cell_frequencies_owen
#'
#' @description
#' # Combining all cells racing lines after binerization (Owen's way)
#'     TODO make the high-points into single points in order to makae the yaxis of the final plot more intuitive.
#'
#' @param list_of_cell_position_frequencies list_of_cell_position_frequencies produced by racing line builder.
#'
#' @return A flat matrix where each column represents a position along a trajectory,
#'     and each value is proportional to the number of the sample's cells which are at that position.
#' @export
flatten_cell_frequencies_owen <- function(list_of_cell_position_frequencies) {
  all_patient_cells_scored_flat <- list()
  length_of_list <- length(list_of_cell_position_frequencies)
  for(i in 1:length_of_list){
    location_of_highpoint <- as.numeric(colSums(list_of_cell_position_frequencies[[i]]) == max(colSums(list_of_cell_position_frequencies[[i]])))
    all_patient_cells_scored_flat[[i]] <- location_of_highpoint
  }
  flat_matrix <- do.call(rbind, all_patient_cells_scored_flat)
  return(flat_matrix)
}
