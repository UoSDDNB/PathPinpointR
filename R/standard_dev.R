#' @title standard_dev
#'
#' @description
#' Calculate the standard deviaiton in the predicted position of cells.
#'
#' @param ppr An object of class ppr.
#'
#' @return the standard deviaiton of the predicted position of cells.
#'
#' @export
#'

calculate_variance <- function(ppr) {
  # Extract the cells_flat matrix from the ppr object
  cells_flat <- ppr$cells_flat
  # Find the pseudo time index for the max value in each cell
  max_values <- apply(cells_flat, 1, which_mid_max)
  #calulate the standard deviation of the pseudo time index
  return(sd(max_values))
}

