#' @title predict_position
#'
#' @description
#' Produces an estimate for the position on trajectory of each cell of a sample.
#' This is later aggregated to estimate samples position along the trajectory.
#'
#' @param sample_sce A Single Cell Experiment object,
#' containing a matrix of your samples binary gene expression,
#' which has been filtered to only include switching genes,
#' using PathPinPointR::subset_switching_genes().
#' @param switching_genes Genes which switch through the trajectory,
#' as identified by GeneSwitches.
#'
#' @return A list of matrices:
#' A matrix for each cell:
#' the columns represent progress through a trajectory,
#' the rows represent genes,
#' values indicate a likely position of the cell upon the trajectory,
#' based that genes bianrized expression.
#' @export
#'
predict_position <- function(sample_sce, switching_genes) {

  # Check if sample_sce contains binarized expression data
  if (is.null(sample_sce@assays@data@listData$binary)) {
    # If binarized expression data is missing, display a message and stop.
    stop("\n  Binarized expression data not found in \"",
         deparse(substitute(sample_sce)),
         "\" \n  You must run `GeneSwitches::binarize_exp()` first.")
  } else {
    reduced_binary_counts_matrix <- sample_sce@assays@data@listData$binary
  }
  ## The final output will be a ppr_obj (list) comprised of 3 objects.
  # Make the list of length 3, and name the objects
  ppr_obj <- vector("list", 4)
  names(ppr_obj) <- c("genomic_expression_traces",
                      "cells_flat",
                      "sample_flat",
                      "sd")
  # Assign the ppr_obj class attribute to the list
  class(ppr_obj) <- "PPR_OBJECT"

  ## Reorder switching_genes
  # because the code relies on the rownames and idicies of the genes,
  # in reduced_binary_counts_matrix and switching_genes matching.
  switching_genes <- switching_genes[rownames(reduced_binary_counts_matrix),]

  ## Input values:
  number_of_cells <- dim(reduced_binary_counts_matrix)[2]
  # Should we get the number of switching genes from here?
  # or dim(reduced_binary_counts_matrix)[1]
  # dim(reduced_binary_counts_matrix)[1] < dim(switching_genes)[1] is possible.
  number_of_switching_genes <- dim(switching_genes)[1]


  # Different datasets store this information in different columns.. watch out.
  # Maybe this should be an input to the function?
  # as.numeric is needed/notneeded depending on dataset.
  switching_time <- as.numeric(switching_genes$switch_at_timeidx)
  switching_direction <- switching_genes$direction

  # Building an empty genomic_expression_traces list.
  # (faster than building it dynamically.)
  all_patients_cells_scored <- vector("list", number_of_cells)
  names(all_patients_cells_scored) <- colnames(reduced_binary_counts_matrix)

  # Loop through all cells,
  # making matrices for each, which represent likely position of a cell (c),
  # on a trajectory based on the expression of each gene.
  for (c in 1:number_of_cells) {
    # Build the matrix of 0's
    # which has genes as rows and pseudotime indecies as columns.
    genomic_expression_mat <- matrix(0,
                                     nrow = number_of_switching_genes, 
                                     ncol = 100)

    rownames(genomic_expression_mat) <- rownames(reduced_binary_counts_matrix)
    bin_exp_c <- reduced_binary_counts_matrix[, c]

    up_indices <- which(bin_exp_c == 1 & switching_direction == "up")

    down_indices <- which(bin_exp_c == 0 & switching_direction == "down")

    not_up_indices <- which(bin_exp_c == 0 & switching_direction == "up")

    not_down_indices <- which(bin_exp_c == 1 & switching_direction == "down")

    for (i in up_indices) {
      genomic_expression_mat[i, switching_time[i]:100] <- 1
    }

    for (i in down_indices) {
      genomic_expression_mat[i, switching_time[i]:100] <- 1
    }

    for (i in not_up_indices) {
      genomic_expression_mat[i, 1:switching_time[i]] <- 1
    }

    for (i in not_down_indices) {
      genomic_expression_mat[i, 1:switching_time[i]] <- 1
    }

    all_patients_cells_scored[[c]] <- genomic_expression_mat
  }

  ppr_obj$genomic_expression_traces <- all_patients_cells_scored

  ### GENOMIC EXPRESSION TRACES CREATED
  # Now flatten:

  # Use lapply to calculate column sums for each matrix
  ppr_obj$cells_flat <- do.call(rbind,
                                lapply(all_patients_cells_scored,
                                       colSums))
  rownames(ppr_obj$cells_flat) <- names(all_patients_cells_scored)
  # Combine each cells column sums into a single flat matrix.
  ppr_obj$sample_flat <- matrix(colSums(ppr_obj$cells_flat), nrow = 1)

  ## Calculate the standard deviation of the predicted position of cells.
  ppr_obj$sd <- sd(apply(ppr_obj$cells_flat, 1, which_mid_max))

  return(ppr_obj)
}
