# R/ml_utils.R

#' @title Prepare Training Data from SingleCellExperiment
#'
#' @description
#' Extracts binary expression matrix (features) and pseudotime (target) from
#' a SingleCellExperiment object for ML model training.
#'
#' @param reference_sce A SingleCellExperiment object containing binarized
#' expression data and pseudotime information.
#' @param gene_features Character vector of gene names to use as features.
#' If NULL, uses all genes in the binary expression matrix.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{X}: Matrix of features (cells × genes)
#'   \item \code{y}: Vector of target values (pseudotime)
#'   \item \code{cell_names}: Character vector of cell names
#' }
#'
#' @keywords internal
prepare_training_data <- function(reference_sce, gene_features = NULL) {
  
  # Check for binarized expression data
  if (is.null(reference_sce@assays@data@listData$binary)) {
    stop("\n  Binarized expression data not found in reference_sce.",
         "\n  You must run `GeneSwitches::binarize_exp()` first.")
  }
  
  # Check for pseudotime
  if (is.null(reference_sce@colData$Pseudotime)) {
    stop("\n  Pseudotime not found in reference_sce.",
         "\n  You must run slingshot or provide pseudotime information.")
  }
  
  # Get binary expression matrix
  binary_matrix <- reference_sce@assays@data@listData$binary
  
  # Convert to sparse matrix for memory efficiency (binary data is mostly zeros)
  # Check if Matrix package is available
  if (requireNamespace("Matrix", quietly = TRUE)) {
    # Convert to sparse matrix if not already
    if (!inherits(binary_matrix, "sparseMatrix")) {
      binary_matrix <- Matrix::Matrix(binary_matrix, sparse = TRUE)
    }
  } else {
    # Fallback to dense matrix if Matrix package not available
    if (!is.matrix(binary_matrix)) {
      binary_matrix <- as.matrix(binary_matrix)
    }
  }
  
  # Subset to specified genes if provided
  if (!is.null(gene_features)) {
    # Check that all requested genes are present
    missing_genes <- setdiff(gene_features, rownames(binary_matrix))
    if (length(missing_genes) > 0) {
      warning(paste("Some requested genes not found:", 
                    paste(missing_genes, collapse = ", ")))
    }
    available_genes <- intersect(gene_features, rownames(binary_matrix))
    if (length(available_genes) == 0) {
      stop("No requested genes found in binary expression matrix.")
    }
    binary_matrix <- binary_matrix[available_genes, , drop = FALSE]
  }
  
  # Transpose to get cells × genes format
  # Matrix::t() works for both sparse and dense matrices
  if (requireNamespace("Matrix", quietly = TRUE) && inherits(binary_matrix, "sparseMatrix")) {
    X <- Matrix::t(binary_matrix)
  } else {
    X <- t(binary_matrix)
  }
  
  # Get pseudotime values
  y <- as.numeric(reference_sce@colData$Pseudotime)
  
  # Get cell names
  cell_names <- colnames(reference_sce)
  
  # Remove any cells with missing pseudotime
  valid_idx <- !is.na(y)
  if (sum(!valid_idx) > 0) {
    warning(paste("Removing", sum(!valid_idx), "cells with missing pseudotime."))
    X <- X[valid_idx, , drop = FALSE]
    y <- y[valid_idx]
    cell_names <- cell_names[valid_idx]
  }
  
  return(list(
    X = X,
    y = y,
    cell_names = cell_names
  ))
}

#' @title Convert Pseudotime to Index Format
#'
#' @description
#' Converts continuous pseudotime values to 1-100 index format used by
#' PathPinpointR, matching the conversion used in accuracy_test().
#'
#' @param pseudotime Numeric vector of pseudotime values.
#' @param pseudotime_range Numeric vector of length 2: c(min, max) of
#' reference pseudotime range. If NULL, calculated from pseudotime.
#'
#' @return Numeric vector of indices (1-100).
#'
#' @keywords internal
convert_pseudotime_to_index <- function(pseudotime, pseudotime_range = NULL) {
  
  if (is.null(pseudotime_range)) {
    pseudotime_range <- c(min(pseudotime, na.rm = TRUE), 
                         max(pseudotime, na.rm = TRUE))
  }
  
  # Calculate step size (matching accuracy_test logic)
  steptime <- (pseudotime_range[2] - pseudotime_range[1]) / 100
  
  # Convert to indices
  indices <- round((pseudotime - pseudotime_range[1]) / steptime)
  
  # Ensure indices are within 1-100 range
  indices <- pmax(1, pmin(100, indices))
  
  return(indices)
}

#' @title Create PPR_OBJECT from ML Predictions
#'
#' @description
#' Converts ML model predictions (continuous pseudotime) into PPR_OBJECT format
#' compatible with existing PathPinpointR plotting and evaluation functions.
#'
#' @param predicted_pseudotime Numeric vector of predicted pseudotime values.
#' @param cell_names Character vector of cell names (must match length of
#' predicted_pseudotime).
#' @param pseudotime_range Numeric vector of length 2: c(min, max) of
#' reference pseudotime range used for training.
#'
#' @return A PPR_OBJECT containing:
#' \itemize{
#'   \item \code{cells_flat}: Matrix (cells × 100) with predicted positions
#'   \item \code{sample_flat}: Matrix (1 × 100) with aggregated predictions
#' }
#'
#' @keywords internal
create_ppr_from_predictions <- function(predicted_pseudotime, 
                                        cell_names,
                                        pseudotime_range) {
  
  # Convert predictions to indices
  predicted_indices <- convert_pseudotime_to_index(predicted_pseudotime, 
                                                   pseudotime_range)
  
  # Create cells_flat matrix
  n_cells <- length(predicted_indices)
  cells_flat <- matrix(0, nrow = n_cells, ncol = 100)
  rownames(cells_flat) <- cell_names
  
  # Set predicted position to 1 (matching current method's approach)
  for (i in 1:n_cells) {
    idx <- predicted_indices[i]
    cells_flat[i, idx] <- 1
  }
  
  # Create sample_flat (aggregated predictions)
  sample_flat <- matrix(colSums(cells_flat), nrow = 1)
  
  # Create PPR_OBJECT
  ppr_obj <- list(
    cells_flat = cells_flat,
    sample_flat = sample_flat
  )
  class(ppr_obj) <- "PPR_OBJECT"
  
  return(ppr_obj)
}
