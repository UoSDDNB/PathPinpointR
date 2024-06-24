#' accuracy_test
#'
#' @description
#' Use this with "sample" cells taken from the reference.
#' Checks how close the predicted position of the cells is to the true position.
#'
#' @param sample_ppr A "PPR_OBJECT" which originates from the reference,
#' Or a subet of the reference.
#' @param reference_sce A SingleCellExperiment object,
#' containing Pseudotime, binarized and found switching genes.
#' @param plot Logical value. Default is FALSE.
#' If TRUE, a histogram of inaccuracy will be plotted.
#'
#' @return If plot = FALSE, returns a data frame containing accuracy results.
#' If plot = TRUE, returns a the accuracy data frame.
#' @importFrom graphics abline hist text
#' @export
accuracy_test <- function(sample_ppr, reference_sce, plot = TRUE) {
  # Create a data frame to store the accuracy results
  accuracy <- data.frame(
    # Cell names from reference_sce
    cell_names = colnames(reference_sce),
    # True pseudotime values (from reference_sce)
    true_pseudotime = reference_sce@colData$Pseudotime,
    # True time indices
    true_timeIDX = NA,
    # Predicted time indices
    predicted_timeIDX = NA,
    # Placeholder for inaccuracy values
    inaccuracy = NA
  )

  #
  ref_ptime <- reference_sce@colData$Pseudotime

  # Calculate the time step based on the pseudotime range
  # TODO clarify the deinfition of "time step"
  steptime <- (max(ref_ptime) - min(ref_ptime)) / 100

  # Calculate the true time indices based on the pseudotime values
  accuracy$true_timeIDX <- round((ref_ptime -
                                    min(ref_ptime)) / steptime)

  # Predict the time indices for the cells using the sample_ppr
  # Note:
  # This assumes that the names of cells in sample_ppr$cells_flat match the,
  # cell names in reference_sce
  accuracy$predicted_timeIDX[match(names(apply(sample_ppr$cells_flat, 1, which_mid_max)), accuracy$cell_names)] <- apply(sample_ppr$cells_flat, 1, which_mid_max) 
  # Faster method which only works with reference being used as sample.
  #accuracy$predicted_timeIDX <- max.col(sample_ppr$cells_flat, "first")

  # Calculate the accuracy
  # (the absolute difference between the true and predicted time indices)
  accuracy$inaccuracy <- abs(accuracy$true_timeIDX - accuracy$predicted_timeIDX)

  if (plot) {
    invisible({
      hist_plot <- hist(accuracy$inaccuracy,
                        breaks = seq(0, max(accuracy$inaccuracy), by = 1),
                        main = "Histogram of Inaccuracy",
                        xlab = "Inaccuracy")
      mean_inaccuracy <- mean(accuracy$inaccuracy, na.rm = TRUE)
      abline(v = mean_inaccuracy, col = "red", lwd = 1)
      text(mean_inaccuracy,
           max(hist_plot$counts) * 0.9,
           labels = paste("Mean =",
                          round(mean_inaccuracy, 2)),
           adj = c(0.5, 0),
           col = "red")
      hist_plot
    })
  } else {
    return(accuracy)
  }
}
