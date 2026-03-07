#' pt_bias
#'
#' @description
#' This function plots the accuracy of PPR predictions across pseudo-time. 
#' 
#' @param ppr A "PPR_OBJECT" which originates from the reference,
#' Or a subset of the reference.
#' @param reference_sce A SingleCellExperiment object.
#' This object should contain:
#' - `Pseudotime`: A numeric vector representing the progression of cells along a trajectory.
#' - `Binarized`: A matrix or data frame indicating the binary expression status of genes.
#' @return A ggplot object showing a scatter plot of predicted vs. true pseudotime indices, 
#' with a reference line indicating perfect prediction.
#' These elements should originate from the reference data.
#'
#' @return returns a ggplot object showing the bias of PPR predictions across pseudo-time.
#' @importFrom ggplot2 ggplot aes geom_point labs geom_abline xlim ylim theme
#' @importFrom ggplot2 ggplot aes geom_point labs geom_abline xlim ylim theme margin
#' @export
#' 
pt_bias <- function(ppr, reference_sce) {
  # Get the accuracy data
  accuracy <- accuracy_test(ppr, reference_sce, plot = FALSE)
  
  # Create a ggplot object to show the bias of PPR predictions across pseudo-time
  p <- ggplot(accuracy, aes(x = true_timeIDX, y = predicted_timeIDX)) +
              geom_point() +
              labs(x = "True Pseudotime IDX", y = "Predicted Pseudotime IDX")  +
              geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
              xlim(0, NA) +
              ylim(0, NA) +
              theme(plot.margin = margin(l = 20, b = 20, r = 20, t = 20))
  
  return(p)
}