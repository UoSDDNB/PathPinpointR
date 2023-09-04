#' Title flatten_cell_frequencies_moi
#'
#' @description
#' Combining all cells racing lines with out binerization (Moi's way).
#'
#'
#' @param list_of_cell_position_frequencies list_of_cell_position_frequencies produced by racing line builder.
#'
#' @return A flat matrix where each column represents a position along a trajectory,
#'     and each value is proportional to the number of the sample's cells which are at that position.
#' @export
flatten_cell_frequencies_moi <- function(list_of_cell_position_frequencies) {
  #making an empty flat matrix
  all_patient_cells_scored_flat <- matrix(0, nrow = 1, ncol = 100)

  # for every cells matrix calculate the colsums and add them to the flat matrix
  for (i in 1:length(list_of_cell_position_frequencies)) {
    all_patient_cells_scored_flat <- all_patient_cells_scored_flat + colSums(list_of_cell_position_frequencies[[i]])
  }
  return(all_patient_cells_scored_flat)
}




