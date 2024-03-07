#' @title Identify and Visualise each cells position.
#'
#' @description
#' Produces a plot for each cell which helps visualize how PathPinPointR is predicting the cells position.
#'
#' @param reference.sg A selection of switching genes which are evenly distributed through pseudo-time.
#' @param lines logical, do you want to include lines to indicate the predicted position of a cell
#' @param reference_reduced a matrix of your sample's binary gene expression.
#' @param cell_id The index or name of the cell of interest
#'
#' @return Timeline plot of selected cell
#' @import ggplot2
#' @import ggrepel
#'
#'
#' @export
#'

ppr_timeline_plot <- function(reference.sg, lines = TRUE, reference_reduced = NULL, cell_id = 1) {

  # Convert reference.sg to a data frame
  reference.sg <- as.data.frame(reference.sg)

  # Add a new column direction_num and set it to -1
  reference.sg$direction_num <- -1

  # If "up" is present in the direction column, set direction_num to 1
  if ("up" %in% reference.sg$direction) {
    reference.sg[reference.sg$direction == "up", ]$direction_num <- 1
  }

  # Generate pseudotime data based on the specified parameters
    timeidx_df <- data.frame(timeidx_range = c(0, 25, 50, 75, 100), timeidx_labs = c(0, 25, 50, 75, 100))


  # Create the initial ggplot object with x and y aesthetics, color, and labels
  timeline_plot <- ggplot(reference.sg, aes(x = switch_at_timeidx, y = pseudoR2s * direction_num, label = rownames(reference.sg))) +
    geom_point(size = 1) + xlab("Time-Index") + ylab("Quality of fitting (R^2)")

  # Add the classic theme to the plot
  timeline_plot <- timeline_plot + theme_classic()

  # Add a horizontal black line for the timeline
  timeline_plot <- timeline_plot + geom_hline(yintercept = 0, color = "black", linewidth = 0.6)

  # Add labels for pseudotime on the plot
  timeline_plot <- timeline_plot + geom_label(data = timeidx_df, aes(x = timeidx_range, y = 0, label = timeidx_labs),
                                    size = (3), color = "black")

  # Add text labels with repulsion to avoid overlap
  timeline_plot <- timeline_plot + geom_text_repel(aes(x = switch_at_timeidx, y = pseudoR2s * direction_num, label = rownames(reference.sg)),
                                         size = 3, show.legend = FALSE)

  # Customize the theme and legend appearance
  timeline_plot <- timeline_plot + theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(10, "pt"),
                               text = element_text(size = 12))


if (lines) {

  if (is.character(cell_id)) {
    cell_idx <- which(colnames(reference_reduced) == cell_id)
  } else if (is.numeric(cell_id)) {
    # Convert numeric to integer if needed
    cell_idx <- as.integer(cell_id)
  } else {
    stop("cell_id must be character or numeric.")
  }

  ## Reorder reference.sg
  # as the code relies on the rownames and idicies of the genes in reference_reduced and reference.sg matching.
  reference_reduced <- reference_reduced[rownames(reference.sg),]

  # Check if row names in reference.sg match those in reference_reduced
  if (!setequal(rownames(reference.sg), rownames(reference_reduced))) {
    stop("Row names in reference.sg do not match those in reference_reduced")
  }

  #loop through all of the genes in reference.sg.
    for (g in 1:dim(reference.sg)[1]) {
      # if G is NOT expressed in  C, and the switch is UP  , then draw the line to the right.
      if ((reference_reduced[rownames(reference.sg)[g], cell_idx] == 0) && (reference.sg$direction[g] == "up")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = 0, xend = reference.sg$switch_at_timeidx[g], y = reference.sg$pseudoR2s[g], yend = reference.sg$pseudoR2s[g]),
                                                      color = "blue", linewidth = 0.6)
      }
      # if G is expressed in      c, and the switch is UP  , then draw the line to the right.
      if ((reference_reduced[rownames(reference.sg)[g], cell_idx] == 1) && (reference.sg$direction[g] == "up")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = reference.sg$switch_at_timeidx[g], xend = 100, y = reference.sg$pseudoR2s[g], yend = reference.sg$pseudoR2s[g]),
                                                      color = "blue", linewidth = 0.6)
      }
      # if G is NOT expressed in  C, and the switch is Down, then draw the line to the Left.
      if ((reference_reduced[rownames(reference.sg)[g], cell_idx] == 0) && (reference.sg$direction[g] == "down")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = reference.sg$switch_at_timeidx[g], xend = 100, y = -reference.sg$pseudoR2s[g], yend = -reference.sg$pseudoR2s[g]),
                                                      color = "blue", linewidth = 0.6)
      }
      # if G is expressed in      C, and the switch is Down, then draw the line to the Left.
      if ((reference_reduced[rownames(reference.sg)[g], cell_idx] == 1) && (reference.sg$direction[g] == "down")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = 0, xend = reference.sg$switch_at_timeidx[g], y = -reference.sg$pseudoR2s[g], yend = -reference.sg$pseudoR2s[g]),
                                                      color = "blue", linewidth = 0.6)
      }
    }

# Title the plot with the name of the cell
timeline_plot <- timeline_plot + ggtitle(colnames(reference_reduced)[cell_idx])


}

# Return the final plot
return(timeline_plot)

}

