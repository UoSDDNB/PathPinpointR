#' @title Identify and Visualise each cells position.
#'
#' @description
#' Produces a plot for each cell which helps visualize how GSS is predicting the cells position.
#'
#' @param reduced_binary_counts_matrix a matrix of your samples binary gene expression.
#' @param gss_genes A selection of switching genes which are evenly distributed through pseudo-time.
#' @param cell_idx The index (should get changed to name) of the cell of interest
#' @param full_time_IDX Do you want the scale to go from Min to Max or from 0-100.
#'
#' @return Timeline plot of selected cell
#' @import ggplot2
#' @import ggrepel
#'
#'
#' @export
#'

racinglines_timeline <- function(gss_genes, reduced_binary_counts_matrix, cell_idx = 1) {

  # Convert gss_genes to a data frame
  gss_genes <- as.data.frame(gss_genes)

  # Check if row names in gss_genes match those in reduced_binary_counts_matrix
  if (!identical(rownames(gss_genes), rownames(reduced_binary_counts_matrix))) {
    stop("Row names in gss_genes do not match those in reduced_binary_counts_matrix")
  }

  ## Reorder gss_genes
  # as the code relies on the rownames and idicies of the genes in reduced_binary_counts_matrix and gss_genes matching.
  gss_genes <- gss_genes[rownames(reduced_binary_counts_matrix),]

  # Add a new column direction_num and set it to -1
  gss_genes$direction_num <- -1

  # If "up" is present in the direction column, set direction_num to 1
  if ("up" %in% gss_genes$direction) {
    gss_genes[gss_genes$direction == "up", ]$direction_num <- 1
  }

  # Generate pseudotime data based on the specified parameters
    timeidx_df <- data.frame(timeidx_range = c(0, 25, 50, 75, 100), timeidx_labs = c(0, 25, 50, 75, 100))


  # Create the initial ggplot object with x and y aesthetics, color, and labels
  timeline_plot <- ggplot(gss_genes, aes(x = switch_at_timeidx, y = pseudoR2s * direction_num, label = rownames(gss_genes))) +
    geom_point(size = 1) + xlab("Time-Index") + ylab("Quality of fitting (R^2)")

  # Add the classic theme to the plot
  timeline_plot <- timeline_plot + theme_classic()

  # Add a horizontal black line for the timeline
  timeline_plot <- timeline_plot + geom_hline(yintercept = 0, color = "black", linewidth = 0.6)

  # Add labels for pseudotime on the plot
  timeline_plot <- timeline_plot + geom_label(data = timeidx_df, aes(x = timeidx_range, y = 0, label = timeidx_labs),
                                    size = (3), color = "black")

  # Add text labels with repulsion to avoid overlap
  timeline_plot <- timeline_plot + geom_text_repel(aes(x = switch_at_timeidx, y = pseudoR2s * direction_num, label = rownames(gss_genes)),
                                         size = 3, show.legend = FALSE)

  # Customize the theme and legend appearance
  timeline_plot <- timeline_plot + theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(10, "pt"),
                               text = element_text(size = 12))




  #loop through all of the genes in gss_genes.
    for (g in 1:dim(gss_genes)[1]) {
      # if G is NOT expressed in  C, and the switch is UP  , then draw the line to the right.
      if ((reduced_binary_counts_matrix[rownames(gss_genes)[g], cell_idx] == 0) && (gss_genes$direction[g] == "up")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = 0, xend = gss_genes$switch_at_timeidx[g], y = gss_genes$pseudoR2s[g], yend = gss_genes$pseudoR2s[g]),
                                                      color = "blue", size = 0.8)
      }
      # if G is expressed in      c, and the switch is UP  , then draw the line to the right.
      if ((reduced_binary_counts_matrix[rownames(gss_genes)[g], cell_idx] == 1) && (gss_genes$direction[g] == "up")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = gss_genes$switch_at_timeidx[g], xend = 100, y = gss_genes$pseudoR2s[g], yend = gss_genes$pseudoR2s[g]),
                                                      color = "blue", size = 0.8)
      }
      # if G is NOT expressed in  C, and the switch is Down, then draw the line to the Left.
      if ((reduced_binary_counts_matrix[rownames(gss_genes)[g], cell_idx] == 0) && (gss_genes$direction[g] == "down")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = gss_genes$switch_at_timeidx[g], xend = 100, y = -gss_genes$pseudoR2s[g], yend = -gss_genes$pseudoR2s[g]),
                                                      color = "blue", size = 0.8)
      }
      # if G is expressed in      C, and the switch is Down, then draw the line to the Left.
      if ((reduced_binary_counts_matrix[rownames(gss_genes)[g], cell_idx] == 1) && (gss_genes$direction[g] == "down")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = 0, xend = gss_genes$switch_at_timeidx[g], y = -gss_genes$pseudoR2s[g], yend = -gss_genes$pseudoR2s[g]),
                                                      color = "blue", size = 0.8)
      }
    }

# Title the plot with the name of the cell
timeline_plot <- timeline_plot + ggtitle(colnames(reduced_binary_counts_matrix)[cell_idx])

# Return the final plot
return(timeline_plot)

}

