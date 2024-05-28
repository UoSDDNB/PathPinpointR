#' @title Identify and Visualise each cells position.
#'
#' @description Produces a plot for a given cell,
#' this helps visualize predicted position of the selected cell.
#'
#' @param switching_genes A selection of switching genes.
#' @param genomic_expression_traces Logical,
#' Do you want lines which indicate the predicted position of the selected cell.
#' @param reduced_sce a SingleCellExperiment object,
#' that you have already run PathPinPointR::subset_switching_genes() on.
#' @param cell_id The index or name of the cell of interest
#'
#' @return Timeline plot of selected cell
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme_classic geom_hline
#' geom_label element_blank unit element_text aes_string geom_segment ggtitle
#' @importFrom ggrepel geom_text_repel
#'
#'
#' @export
ppr_timeline_plot <- function(switching_genes,
                              genomic_expression_traces = FALSE,
                              reduced_sce = NULL,
                              cell_id = 1) {

  # Convert switching_genes to a data frame
  switching_genes <- as.data.frame(switching_genes)

  # Add a new column direction_num and set it to -1
  switching_genes$direction_num <- -1

  # If "up" is present in the direction column, set direction_num to 1
  if ("up" %in% switching_genes$direction) {
    switching_genes[switching_genes$direction == "up", ]$direction_num <- 1
  }

  # Generate pseudotime data based on the specified parameters
  timeidx_df <- data.frame(timeidx_range = c(0, 25, 50, 75, 100), 
                           timeidx_labs = c(0, 25, 50, 75, 100))


  # Create the initial ggplot object with x and y aesthetics, color, and labels
  timeline_plot <- ggplot(switching_genes,
                          aes(x = switching_genes$switch_at_timeidx,
                              y = switching_genes$pseudoR2s *
                                switching_genes$direction_num,
                              label = rownames(switching_genes))) +
    geom_point(size = 1) +
    xlab("Pseudotime-Index") +
    ylab("Quality of fitting (R^2)")

  # Add the classic theme to the plot
  timeline_plot <- timeline_plot + theme_classic()

  # Add a horizontal black line for the timeline
  timeline_plot <- timeline_plot +
    geom_hline(yintercept = 0,
               color = "black",
               linewidth = 0.6)

  # Add labels for pseudotime on the plot
  timeline_plot <- timeline_plot +
    geom_label(data = timeidx_df,
               aes(x = timeidx_df$timeidx_range,
                   y = 0,
                   label = timeidx_df$timeidx_labs),
               size = (3),
               color = "black")

  # Add text labels with repulsion to avoid overlap
  timeline_plot <- timeline_plot +
    geom_text_repel(aes(x = switch_at_timeidx,
                        y = pseudoR2s * direction_num,
                        label = rownames(switching_genes)),
                    size = 3, 
                    show.legend = FALSE)

  # Customize the theme and legend appearance
  timeline_plot <- timeline_plot +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.key.size = unit(10, "pt"),
          text = element_text(size = 12))


  if (genomic_expression_traces) {

    if (is.character(cell_id)) {
      cell_idx <- which(colnames(reference_reduced) == cell_id)
    } else if (is.numeric(cell_id)) {
      # Convert numeric to integer if needed
      cell_idx <- as.integer(cell_id)
    } else {
      stop("cell_id must be character or numeric.")
    }

    # Check if reduced_sce contains binarized expression data
    if (is.null(reduced_sce@assays@data@listData$binary)) {
      # If binarized expression data is missing, display a message and stop.
      stop("\n  Binarized expression data not found in \"",
           deparse(substitute(sce)),
           "\" \n  You must run `GeneSwitches::binarize_exp()` first.")
    }
  # Extract the binary counts matrix from the reduced_sce object
  bin_matrix <- reduced_sce@assays@data@listData$binary

    # Check if row names in switching_genes match bin_matrix
    if (!setequal(rownames(switching_genes),
                   rownames(bin_matrix))) {
      stop("Row names in switching_genes do not match those in bin_matrix")
    }

    ## Reorder switching_genes
    # as the code relies on the rownames and idicies of the genes in bin_matrix,
    # and switching_genes matching.
    switching_genes <- switching_genes[rownames(bin_matrix),]

    # loop through all of the genes in switching_genes.
    for (g in 1:dim(switching_genes)[1]) {

      # IF G is NOT expressed in C,
      # AND the switch is UP,
      # then draw the line to the right.
      if ((bin_matrix[rownames(switching_genes)[g], cell_idx] == 0) &&
            (switching_genes$direction[g] == "up")) {
        timeline_plot <- timeline_plot +
          geom_segment(aes_string(x = 0,
                       xend = switching_genes$switch_at_timeidx[g],
                       y = switching_genes$pseudoR2s[g],
                       yend = switching_genes$pseudoR2s[g]),
                       color = "blue", linewidth = 0.6)
      }
      # IF G is expressed in C,
      # AND the switch is UP, 
      # then draw the line to the right.
      if ((bin_matrix[rownames(switching_genes)[g], cell_idx] == 1) && (switching_genes$direction[g] == "up")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = switching_genes$switch_at_timeidx[g], xend = 100, y = switching_genes$pseudoR2s[g], yend = switching_genes$pseudoR2s[g]),
                                                      color = "blue", linewidth = 0.6)
      }
      # if G is NOT expressed in  C, and the switch is Down, then draw the line to the Left.
      if ((bin_matrix[rownames(switching_genes)[g], cell_idx] == 0) && (switching_genes$direction[g] == "down")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = switching_genes$switch_at_timeidx[g], xend = 100, y = -switching_genes$pseudoR2s[g], yend = -switching_genes$pseudoR2s[g]),
                                                      color = "blue", linewidth = 0.6)
      }
      # if G is expressed in      C, and the switch is Down, then draw the line to the Left.
      if ((bin_matrix[rownames(switching_genes)[g], cell_idx] == 1) && (switching_genes$direction[g] == "down")) {
        timeline_plot <- timeline_plot + geom_segment(aes_string(x = 0, xend = switching_genes$switch_at_timeidx[g], y = -switching_genes$pseudoR2s[g], yend = -switching_genes$pseudoR2s[g]),
                                                      color = "blue", linewidth = 0.6)
      }
    }

# Title the plot with the name of the cell
timeline_plot <- timeline_plot + ggtitle(colnames(bin_matrix)[cell_idx])

}

# Return the final plot
return(timeline_plot)

}

