#' @title Identify the "racing lines"
#'
#' @description
#' Produces an estimate for the position on trajectory of each gene in each cell of a sample.
#'     This can be later aggregated to estimate the position of the sample along the trajectory.
#'
#' @param reduced_binary_counts_matrix a matrix of your samples binary gene expression.
#' @param reference.sg Genes which switch through the trajectory as identified by GeneSwitches.
#'
#' @return A list of matrices: A matrix for each cell where the columns represent progress through a trajectory,
#'     and the rows represent genes, values indicate a likely position of the cell upon the trajectory based that genes bianrized expression.
#' @export
#'
pppr_predict_position <- function(reduced_binary_counts_matrix,reference.sg) {

  ## The final output will be a pppr_obj (list) comprised of 3 objects.
  # Make the list of length 3, and name the objects
  pppr_obj <- vector("list", 3)
  names(pppr_obj) <- c("racing_lines","cells_flat","sample_flat")
  # Assign the pppr_obj class attribute to the list
  class(pppr_obj) <- "PPPR_OBJECT"

  ## Reorder reference.sg
  # as the code relies on the rownames and idicies of the genes in reduced_binary_counts_matrix and reference.sg matching.
  reference.sg <- reference.sg[rownames(reduced_binary_counts_matrix),]

  ## Input values:
  number_of_cells <- dim(reduced_binary_counts_matrix)[2]
  # Should we get the number of switching genes from here? or dim(reduced_binary_counts_matrix)[1]
  # It is possible for dim(reduced_binary_counts_matrix)[1] <  dim(reference.sg)[1]
  number_of_switching_genes <- dim(reference.sg)[1]


  # Different datasets store this information in different columns.. watch out.
  # Maybe this should be an input to the function?
  # as.numeric is needed/notneeded depending on dataset.
  switching_time <- as.numeric(reference.sg$switch_at_timeidx)
  switching_direction <- reference.sg$direction

  # Building the racing lines list. (faster than building it dynamically.)
  all_patients_cells_scored <- vector("list", number_of_cells)
  names(all_patients_cells_scored) <- colnames(reduced_binary_counts_matrix)

  # Loop through all cells making matrices for each,
  #which represent likely position of a cell on a trajectory based on the expression of each gene.
  for (c in 1:number_of_cells) {
    # Build the matrix of 0's which has genes as rows and pseudotime indecies as columns.
   racing_mat <- matrix(0, nrow = number_of_switching_genes, ncol = 100)
   rownames(racing_mat) <- rownames(reduced_binary_counts_matrix)
    binarized_gene_expression_for_cell_c <- reduced_binary_counts_matrix[, c]

    up_indices <- which(binarized_gene_expression_for_cell_c == 1 & switching_direction == "up")
    down_indices <- which(binarized_gene_expression_for_cell_c == 0 & switching_direction == "down")
    not_up_indices <- which(binarized_gene_expression_for_cell_c == 0 & switching_direction == "up")
    not_down_indices <- which(binarized_gene_expression_for_cell_c == 1 & switching_direction == "down")

    for (i in up_indices) {
      racing_mat[i, switching_time[i]:100] <- 1
    }

    for (i in down_indices) {
      racing_mat[i, switching_time[i]:100] <- 1
    }

    for (i in not_up_indices) {
      racing_mat[i, 1:switching_time[i]] <- 1
    }

    for (i in not_down_indices) {
      racing_mat[i, 1:switching_time[i]] <- 1
    }

    all_patients_cells_scored[[c]] <- racing_mat
  }

  pppr_obj$racing_lines <- all_patients_cells_scored

### RACING LINES CREATED
# Now flatten:

# Use lapply to calculate column sums for each matrix
pppr_obj$cells_flat <- do.call(rbind, lapply(all_patients_cells_scored, colSums))
rownames(pppr_obj$cells_flat) <- names(all_patients_cells_scored)
# Combine each cells column sums into a single flat matrix.
pppr_obj$sample_flat <- matrix(colSums(pppr_obj$cells_flat), nrow = 1)

return(pppr_obj)
}
