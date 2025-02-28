#' @title predict_position
#'
#' @description
#' Produces an estimate for the position on trajectory of each cell of a sample.
#' This is later aggregated to estimate samples position along the trajectory.
#'
#' @param sample_sce A Single Cell Experiment object,
#' containing a matrix of your samples binary gene expression.
#' @param switching_genes Genes which switch through the trajectory,
#' as identified by GeneSwitches.
#'
#' @return A list of matrices:
#' A matrix for each cell:
#' the columns represent progress through a trajectory,
#' the rows represent genes,
#' values indicate a likely position of the cell upon the trajectory,
#' based on that gene's binarized expression.
#' 
#' @importFrom stats sd
#' 
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

  ## Create final output
  # TODO: change this to only be cells_flat.
  ppr_obj <- vector("list", 3)
  names(ppr_obj) <- c("genomic_expression_traces",
                      "cells_flat",
                      "sample_flat")
  # Assign the ppr_obj class attribute to the list
  class(ppr_obj) <- "PPR_OBJECT"

  ## Reorder switching_genes
  # because the code relies on the rownames and idicies of the genes in,
  # reduced_binary_counts_matrix and switching_genes matching.
  # TODO: inlcude a check that reduce_counts_matrix() has been run.
  # Something like setequal(rownames(reduced_binary_counts_matrix), rownames(switching_genes))
  switching_genes <- switching_genes[rownames(reduced_binary_counts_matrix),]

  ## Input values:
  number_of_cells <- dim(reduced_binary_counts_matrix)[2]
  number_of_switching_genes <- dim(switching_genes)[1]
  switching_time <- as.numeric(switching_genes$switch_at_timeidx)
  switching_direction <- switching_genes$direction

  # pre allocate cells_flat
  ppr_obj$cells_flat <- matrix(0, nrow = number_of_cells, ncol = 100)
  rownames(ppr_obj$cells_flat) <- colnames(reduced_binary_counts_matrix)


  # build an empty genomic_expression_mat outside of loop for speed.
  # which has genes as rows and pseudotime indecies as columns.
  empty_genomic_expression_mat <- matrix(0,
                                         nrow = number_of_switching_genes, 
                                         ncol = 100)
  # Set the rownames of the matrix to the gene names.
  rownames(empty_genomic_expression_mat) <- rownames(reduced_binary_counts_matrix)

  # Print a message to the console
  cat("Predicting position of cells... \n")

  # Loop through all cells,
  # Making a matrix for each cell, 
  # the genomic_expression_mat represents the likely position of the cell on the trajectory.
  # The colSums of the matrix is then used to create a row in the cells_flat matrix.
  # cells_flat has cells as rows and pseudotime indecies as columns.
  # the colSums of cells_flat represents the predicted position of the sample. (sample_flat)
  for (c in 1:number_of_cells) {
    # no need to rest the matrix to 0's, as it is done outside the loop.
    genomic_expression_mat <- empty_genomic_expression_mat

    # extract the binary expression of the cell
    bin_exp_c <- reduced_binary_counts_matrix[, c]

    ## find the indices of the genes where the cell has passes their switching point.
    # genes where they are expressed in the cell and switch up during the trajectory.
    up_indices <- which(bin_exp_c == 1 & switching_direction == "up")
    # genes where they are not expressed in the cell and switch down during the trajectory.
    down_indices <- which(bin_exp_c == 0 & switching_direction == "down")

    ## find the indices of the genes where the cell has Not yet passed their switching point.
    # genes where they are not expressed in the cell and switch up during the trajectory.
    not_up_indices <- which(bin_exp_c == 0 & switching_direction == "up")
    # genes where they are expressed in the cell and switch down during the trajectory.
    not_down_indices <- which(bin_exp_c == 1 & switching_direction == "down")

    # set the values of the matrix to 1
    # where the cell has passed the switching point, for the appropriate genes.
    for (i in up_indices) {
      genomic_expression_mat[i, switching_time[i]:100] <- 1
    }
    for (i in down_indices) {
      genomic_expression_mat[i, switching_time[i]:100] <- 1
    }

    # set the values of the matrix to 1
    # where the cell has Not yet passed the switching point, for the appropriate genes.
    for (i in not_up_indices) {
      genomic_expression_mat[i, 1:switching_time[i]] <- 1
    }

    for (i in not_down_indices) {
      genomic_expression_mat[i, 1:switching_time[i]] <- 1
    }

    # assign the colSums of genomic_expression_mat to ppr_obj$cells_flat
    ppr_obj$cells_flat[c, ] <- colSums(genomic_expression_mat)

    # print a message to the console
    if (c %% 100 == 0) {
      cat("\rPredicted ", c, "/", number_of_cells)
       flush.console()  # Ensures output updates in buffered consoles
    }
  }

  # print a message to the console
  cat("\rPredicted ", c, "/", number_of_cells, "\n")

  # dont assign the genomic_expression_traces to the ppr_obj
  # doing so uses too much RAM.
  # TODO: this should be removed from PPR_OBJECT, and this funcion.
  ppr_obj$genomic_expression_traces <- NULL

  # Combine each cells column sums into a single flat matrix.
  # TODO: this can also be removed from the PPR_OBJECT, and this funcion.
  ppr_obj$sample_flat <- matrix(colSums(ppr_obj$cells_flat), nrow = 1)

  ## Calculate the standard deviation of the predicted position of cells.
  #ppr_obj$sd <- sd(apply(ppr_obj$cells_flat, 1, which_mid_max))
  # commented for computational efficiency

  # print a message to the console
  cat("Prediction complete \n") 
  return(ppr_obj)
}