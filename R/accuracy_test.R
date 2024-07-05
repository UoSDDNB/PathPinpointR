#' accuracy_test
#'
#' @description
#' Use this with "sample" cells taken from the reference.
#' Checks how close the predicted position of the cells is to the true position.
#'
#' @param ppr A "PPR_OBJECT" which originates from the reference,
#' Or a subet of the reference.
#' @param reference_sce A SingleCellExperiment object,
#' containing Pseudotime, binarized and found switching genes,
#' originating from the reference data.
#' @param plot Logical value. Default is FALSE.
#' If TRUE, a histogram of inaccuracy will be plotted.
#' @param random Logical value. Default is TRUE.
#' If TRUE, a histogram of random inaccuracy will be plotted. 
#' random_inncauracy is the difference of true_time_idx to a random prediction.
#'
#' @return If plot = FALSE, returns a data frame containing accuracy results.
#' If plot = TRUE, returns a the accuracy data frame.
#' @importFrom graphics abline hist text
#' @export
accuracy_test <- function(ppr, reference_sce, plot = TRUE, random = FALSE) {

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
    inaccuracy = NA,
    # Place holder for random prediction
    random_prediction = NA,
    # Placeholder for random inncauracy
    random_inaccuracy = NA
  )

  # Get the pseudotime values from the reference
  ref_ptime <- reference_sce@colData$Pseudotime

  # Calculate the time step based on the pseudotime range
  # TODO clarify the deinfition of "time step"
  steptime <- (max(ref_ptime) - min(ref_ptime)) / 100

  # Calculate the true time indices based on the pseudotime values
  accuracy$true_timeIDX <- round((ref_ptime - min(ref_ptime)) / steptime)

  # Predict the time indices for the cells using the ppr
  # Note:
  # This assumes that the names of cells in ppr$cells_flat match the,
  # cell names in reference_sce
  accuracy$predicted_timeIDX[match(names(apply(ppr$cells_flat,
                                               1,
                                               which_mid_max)),
                                   accuracy$cell_names)] <-
    apply(ppr$cells_flat,
          1,
          which_mid_max)

  # Calculate the accuracy
  # (the absolute difference between the true and predicted time indices)
  accuracy$inaccuracy <- abs(accuracy$true_timeIDX - accuracy$predicted_timeIDX)

  if (random) {
    ### Calculate the accuracy of a random prediction
    # produce random prediction values
    accuracy$random_prediction <- c(sample(0:100,
                                           nrow(ppr$cells_flat),
                                           replace = TRUE),
      rep(NA, nrow(accuracy) - nrow(ppr$cells_flat))
    )


    # calculate the inncauracy of the random prediction
    accuracy$random_inaccuracy <-
      abs(accuracy$true_timeIDX - accuracy$random_prediction)
  }

  if (plot) {
    invisible({
      hist_plot <- hist(accuracy$inaccuracy,
                        breaks = seq(0, 100, by = 1),
                        col = rgb(1, 0 ,0 , 0.5),
                        main = "Histogram of Accuracy",
                        xlab = "Distance from True Time Index")
      mean_inaccuracy <- mean(accuracy$inaccuracy, na.rm = TRUE)
      abline(v = mean_inaccuracy, col = "black", lwd = 1)
      text(x = mean_inaccuracy,
           max(hist_plot$counts) * 0.9,
           labels = paste("Mean =",
                          round(mean_inaccuracy, 2)),
           adj = c(-0.1, 0),
           col = "black")
      if (random) {
        hist(accuracy$random_inaccuracy,
             col = rgb(0, 0, 1, 0.4),
             breaks = seq(0, 100, by = 1),
             add = TRUE)
        legend_colors <- c(rgb(1, 0, 0, 1), rgb(0, 0, 1, 0.4))
        legend_labels <- c("PathPinpointR predicitons", "Random predictions")
        mean_random_inaccuracy <- mean(accuracy$random_inaccuracy, na.rm = TRUE)
        abline(v = mean_random_inaccuracy,
               col = "black",
               lwd = 1)
        text(x = mean_random_inaccuracy,
             max(hist_plot$counts) * 0.85,
             labels = paste("Mean =",
                            round(mean_random_inaccuracy, 2)),
             adj = c(-0.1, 0),
             col = "black")
        text(x = mean_random_inaccuracy,
             max(hist_plot$counts) * 0.87,
             labels = paste("Sum =",
                            sum(accuracy$random_inaccuracy, na.rm = TRUE)),
             adj = c(-0.1, 0),
             col = "black")
        legend("topright", # Position of the legend
               legend = legend_labels, # Labels
               fill = legend_colors) # Colors
        text(x = mean_inaccuracy,
             max(hist_plot$counts) * 0.92,
             labels = paste("Sum =",
                            sum(accuracy$inaccuracy, na.rm = TRUE)),
             adj = c(-0.1, 0),
             col = "black")
      }
      hist_plot
    })
  } else {
    return(accuracy)
  }
}
