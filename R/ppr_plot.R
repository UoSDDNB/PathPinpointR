#' @title ppr_plot
#'
#' @description Base plot for ppr plotting.
#'
#' @return an empty plotting space to be filled by other ppr plotting functions.
#'
#' @export

ppr_plot <- function() {
  plot(x = 1:100,
       ylim = c(0, 115),
       type = "n", # No points
       xlab = "Pseudo-Time Index",
       ylab = "PPR Score",
       main = "Predicted Positions",
       bty = "n") # Remove box around plot
}

#' @title sample_prediction
#'
#' @description Plots the predicted position of your sample.
#'
#' @param sample_ppr A PPR_OBJECT of the sample you wish to plot
#' @param col The colour that you'd like
#' @param label string that you would like to assign as the label to the line.
#'
#' @return an overlay for the ppr plot,
#' which highlights the probable position of your sample in pseudotime.
#'
#' @importFrom graphics segments text lines
#' @export

sample_prediction <- function(sample_ppr,
                              col = "red",
                              label = "sample name") {

  # check that sample_ppr is a PPR_OBJECT
  if (class(sample_ppr) != "PPR_OBJECT") {
    stop("sample_ppr must be a PPR_OBJECT")
  }

  # extract the sample_flat from the PPR_OBJECT
  sample_flat <- sample_ppr$sample_flat

  # identify the index of the maximum value in the sample
  # aka the point PPR predicts the sample to be in pseudotime
  mid_max_idx <- which_mid_max(colSums(sample_flat))

  # Convert the samples raw ppr values to percentages
  y_vals <- (sample_flat / max(sample_flat) * 100)

  # Plot the ppr score for each pseudo-time index
  lines(x = 1:100, y = y_vals, col = col)

  # Add a dashed line to indicate the maximum
  segments(mid_max_idx,
           -3.9,
           mid_max_idx,
           y_vals[mid_max_idx],
           lwd = 1,
           lty = 2,
           col = col)

  # Add sample name as a label
  text(x = mid_max_idx,
       y = y_vals[mid_max_idx],
       labels = label,
       col = col,
       pos = 3)

}

#' @title reference_idents
#'
#' @description
#' Overlay to ppr plot visualising the pseudotime of reference cells,
#' subset by ident.
#'
#' @param reference_sce A SingleCellExperiemnt object of the ref trajectory
#' @param ident selected column from reference_sce colData.
#'
#' @return overlay for the plot_position function with reference idents.
#'
#' @importFrom SummarizedExperiment colData
#' @export

reference_idents <- function(reference_sce, ident) {

  ## checks

  # check that reference_sce is a SingleCellExperiment object
  if (class(reference_sce) != "SingleCellExperiment") {
    stop("reference_sce must be a SingleCellExperiment object")
  }

  # check that refence_sce has a column named Pseudotime
  if (!"Pseudotime" %in% colnames(colData(reference_sce))) {
    stop("reference_sce must have a column named 'Pseudotime'")
  }

  # check that ident is a column in colData(reference_sce)
  if (!ident %in% colnames(colData(reference_sce))) {
    stop("ident must be a column in colData(reference_sce)")
  }

  ## Create a data frame to store the pseudotime & idents
  ref_lab_df <- data.frame(
    # Cell names from reference_sce
    cell_names = colnames(reference_sce),
    # True pseudotime values (from reference_sce)
    true_pseudotime = colData(reference_sce)[, "Pseudotime"],
    # True time indices
    true_timeIDX = NA,
    # Reference labels
    idents = droplevels(factor(colData(reference_sce)[, ident]))
  )

  ## Calculate the true time indices based on the pseudotime values
  # Get the pseudotime values from the reference
  ref_ptime <- colData(reference_sce)[, "Pseudotime"]
  # divide the pseudotime range into 100 even intervals
  pt_intervals <- (max(ref_ptime) - min(ref_ptime)) / 100
  # Calculate the true time indices based on the pseudotime values
  ref_lab_df$true_timeIDX <- round((ref_ptime - min(ref_ptime)) / pt_intervals)


  ## Plot the reference idents

  # count levels of idents
  n_ident_levels <- length(levels(ref_lab_df$idents))

  #check that the number of levels is greater than 0
  if (n_ident_levels < 1) {
    stop("idents has no levels.")
  }

  # check that the number of levels is less than 10
  # if it is greater than 10, continue but produce a warning
  if (n_ident_levels > 9) {
    warning("Too many reference labels, only <10 supported.")
  }

  # add a horizontal line at y=105
  # to separate the boxplot from the rest of the plot
  abline(h = 102.5, col = "black", lty = 1)

  # Define a colour for each level of idents
  box_colours <- gray(rep(0.5, length.out = n_ident_levels))

  # Generate y positions for each boxplot
  y_positions <- seq(to = 103.5,
                     from = 118.5,
                     length.out = n_ident_levels)

  boxplot(ref_lab_df$true_timeIDX ~ ref_lab_df$idents,
          add = TRUE,              # Overlay the boxplot on the existing plot
          horizontal = TRUE,       # Horizontal to align with the x-axis
          at = y_positions,        # Y position for all boxplots
          boxwex = 1,              # Width of the boxes
          col = box_colours,       # Use different colours for each label
          outline = FALSE,         # Remove Outliers
          yaxt = "n",              # supress y-axis label
          frame.plot = FALSE)      # Remove the frame around the plot

  # add the labels for the boxplot
  axis(side = 2,                           # Add labels to the y-axis
       at = y_positions,                    # Position labels according to boxplot positions
       labels = levels(ref_lab_df$idents),  # categorical labels
       las = 1,                             # Rotate the labels to be vertical
       cex.axis = 0.8)                      # Reduce the size of the labels

}


#########

#' @title switching_times
#'
#' @description Adds gene switching points to a PPR plot.
#' The function overlays dashed lines at the pseudotime points where genes
#' switch on or off, along with labels for each gene.
#'
#' @param genes A character vector of gene names you wish to plot.
#' @param switching_genes as output from GeneSwitches::find_switch_logistic_fastglm()
#' 
#' @return Adds dashed lines and gene names to the PPR plot, indicating
#' pseudotime switching points of interest.
#'
#' @importFrom graphics segments text
#' @export
#'

switching_times <- function(genes, switching_genes) {

  # Loop through each gene to add switching points
  for (gene_name in genes) {
    idx <- switching_genes[gene_name, "switch_at_timeidx"]

    # randomise the height of the gene label
    r_num <- sample(2:16, 1) / 2

    # Add a dashed line at the gene's switching index
    segments(idx,
             -4,       # start y-axis
             idx,
             0 + r_num,       # End at -0.5 for a slight offset
             lwd = 1,    # Line width
             lty = 2)    # Dashed line

    # Add gene name label at the switching point
    text(x = idx + 3,             # Slightly offset from the segment
         y = 0.3 + r_num ,                   # Position near the top
         labels = gene_name,       # Label with gene name
         srt = -20,               # Slanted text
         cex = 0.86,              # Font size
         pos = 2)                 # Left-aligned
  }
}



#' @title metrics
#'
#' @description Adds statistical metrics (Standard Deviation, Z-score, p-value) 
#' for a sample to the PPR plot. The metrics are displayed at the sample's 
#' predicted position in pseudotime.
#'
#' @param sample_ppr A PPR_OBJECT containing the PPR data for the sample.
#' Must include the calculated metrics: standard deviation (`sd`), z-score 
#' (`z_score`), and p-value (`p_value`).
#' @param col A character specifying the color of the text.
#' @param sample_flat
#'
#' @return Overlays the PPR plot with statistical information 
#' (SD, Z-score, p-value) for the sample, displayed near the predicted 
#' pseudotime position.
#'
#' @importFrom graphics text
#' @export

metrics <- function(sample_ppr,
                    col = "black",
                    sample_flat) {
  usethis::deprecated("This is a legacy function and will not be maintained actively.")

  # Produce a warning if any metrics are not available
  if (is.null(sample_ppr$sd) ||
        is.null(sample_ppr$z_score) ||
        is.null(sample_ppr$p_value)) {
    warning("One or more metrics are not available for this sample. 
              Please ensure you have run zscore_and_pvalue().")
  }

  # extract the sample_flat from the PPR_OBJECT
  sample_flat <- sample_ppr$sample_flat

  # identify the index of the maximum value in the sample
  # aka the point PPR predicts the sample to be in pseudotime
  mid_max_idx <- which_mid_max(colSums(sample_flat))

  # Y value for the text labels
  # with jitter
  y_val <- sample(12:25, 1)

  # Add Standard Deviation (SD)
  text(x = mid_max_idx,
       y = y_val - 1,     # Position
       labels = paste("SD = ", round(sample_ppr$sd, 2)), # Round
       cex = 0.8,                                        # Font size
       col = col,
       pos = 4)  # Positioned below the point

  # Add Z-score
  text(x = mid_max_idx,
       y = y_val - 2,   # Slightly lower than SD
       labels = paste("Z = ", round(sample_ppr$z_score, 2)),
       cex = 0.8,
       col = col,
       pos = 4)

  # Add p-value
  text(x = mid_max_idx,
       y = y_val - 3,   # Lower than Z-score
       labels = paste("p = ", sample_ppr$p_value),
       cex = 0.8,
       col = col,
       pos = 4)
}

#' @title cell_position_box
#'
#' @description Creates a boxplot overlay showing the predicted pseudotime
#' positions for each cell in a sample.
#'
#' @param sample_ppr A PPR_OBJECT
#' @param col A character vector specifying the colour for the sample.
#' @param points A logical indicating whether to add points to the boxplot.
#'
#' @return Adds boxplots to the PPR plot, visualising the predicted pseudotime
#' positions for each cell in a sample.
#' Points indicate the predicted position of each cell in the sample.
#'
#' @importFrom graphics boxplot stripchart
#' @export
#'

cell_position_box <- function(sample_ppr, col = "gray", points = FALSE) {

  # plot at a random Y within a range:
  y_pos <- sample(25:44, 1)

  # Plot the boxplot
  boxplot(apply(sample_ppr$cells_flat, 1, which_mid_max),
          add = TRUE,               # Overlay the boxplot on the existing plot
          horizontal = TRUE,        # Horizontal to align with the x-axis
          at = y_pos,       
          boxwex = 4,               # Width of the boxes
          col = col, 
          yaxt = "n",               # Suppress y-axis label
          frame.plot = FALSE)       # Remove the frame around the plot

  # Add points to the boxplot
  if (points) {
    stripchart(apply(sample_ppr$cells_flat, 1, which_mid_max),
               method = "jitter",     # Adds jitter to points to avoid overlap
               add = TRUE,            # Add the points to the existing boxplot
               at = y_pos, # Ensure same y-position as the boxplot
               jitter = 0.4,          # Adjust the amount of jitter
               pch = 16,              # Use filled circles for points
               cex = 0.7,             # Reduce the size of the points
               col = "black")             # Use the same colour as the boxplot
  }
}


