#' @title Predict Position Using ML Model
#'
#' @description
#' Predicts the position of cells along a trajectory using a trained machine learning
#' model. Returns results in the same PPR_OBJECT format as the original
#' \code{predict_position} function for compatibility with existing plotting and
#' evaluation functions.
#'
#' @param sample_sce A SingleCellExperiment object containing binarized expression
#' data (from \code{GeneSwitches::binarize_exp()}).
#' @param trained_model A trained model object (from \code{train_ml_model()}).
#' @param verbose Logical. If TRUE (default), print progress messages.
#'
#' @return A PPR_OBJECT containing:
#' \itemize{
#'   \item \code{cells_flat}: Matrix (cells × 100) with predicted positions
#'   \item \code{sample_flat}: Matrix (1 × 100) with aggregated predictions
#' }
#' The output format matches \code{predict_position()} for compatibility with
#' \code{accuracy_test()}, \code{ppr_plot()}, and other PathPinpointR functions.
#'
#' @details
#' This function uses a trained machine learning model (Ridge regression, Random Forest,
#' or XGBoost) to predict continuous pseudotime values from binary gene expression.
#' The predictions are then converted to the 1-100 index format used by PathPinpointR.
#'
#' The function automatically matches genes between the sample and the trained
#' model. If genes are missing, they are set to 0 (not expressed). If extra
#' genes are present, they are ignored.
#'
#' @examples
#' \dontrun{
#' # Train model on reference data
#' trained_model <- train_ml_model(reference_sce, 
#'                                  switching_genes = switching_genes,
#'                                  genes = "switching")
#'
#' # Predict position of sample cells
#' sample_ppr <- predict_position_ml(sample_sce, trained_model)
#'
#' # Evaluate accuracy (if sample is from reference)
#' accuracy_test(sample_ppr, reference_sce, plot = TRUE)
#' }
#'
#' @importFrom glmnet cv.glmnet
#' @export
predict_position_ml <- function(sample_sce, trained_model, verbose = TRUE) {
  
  # Validate trained_model
  if (is.null(trained_model) || !inherits(trained_model, "PPR_ML_MODEL")) {
    stop("trained_model must be a PPR_ML_MODEL object from train_ml_model()")
  }
  
  # Get method from trained model
  method <- trained_model$method
  
  # Check for binarized expression data
  if (is.null(sample_sce@assays@data@listData$binary)) {
    stop("\n  Binarized expression data not found in \"",
         deparse(substitute(sample_sce)),
         "\" \n  You must run `GeneSwitches::binarize_exp()` first.")
  }
  
  # Get binary expression matrix
  binary_matrix <- sample_sce@assays@data@listData$binary
  
  # Convert to sparse matrix for memory efficiency if Matrix package available
  if (requireNamespace("Matrix", quietly = TRUE)) {
    if (!inherits(binary_matrix, "sparseMatrix")) {
      binary_matrix <- Matrix::Matrix(binary_matrix, sparse = TRUE)
    }
  } else {
    if (!is.matrix(binary_matrix)) {
      binary_matrix <- as.matrix(binary_matrix)
    }
  }
  
  # Get required gene features from model
  required_genes <- trained_model$gene_features
  
  # Check which genes are available
  available_genes <- rownames(binary_matrix)
  missing_genes <- setdiff(required_genes, available_genes)
  present_genes <- intersect(required_genes, available_genes)
  
  if (length(present_genes) == 0) {
    stop("None of the required genes are present in sample_sce.")
  }
  
  if (length(missing_genes) > 0 && verbose) {
    warning(paste("Missing", length(missing_genes), 
                  "genes in sample. Setting to 0 (not expressed)."))
  }
  
  # Create feature matrix matching model's gene order
  # Method-specific matrix creation
  n_cells <- ncol(binary_matrix)
  
  if (method == "ridge_regression") {
    # For Ridge: use sparse matrix if available (memory efficient)
    if (requireNamespace("Matrix", quietly = TRUE)) {
      # Create sparse matrix (more memory efficient for binary data)
      X <- Matrix::Matrix(0, nrow = n_cells, ncol = length(required_genes), sparse = TRUE)
      colnames(X) <- required_genes
      rownames(X) <- colnames(binary_matrix)
      
      # Fill in available genes
      if (length(present_genes) > 0) {
        # Extract and transpose subset
        if (inherits(binary_matrix, "sparseMatrix")) {
          X[, present_genes] <- Matrix::t(binary_matrix[present_genes, , drop = FALSE])
        } else {
          X[, present_genes] <- t(binary_matrix[present_genes, , drop = FALSE])
        }
      }
    } else {
      # Fallback to dense matrix
      X <- matrix(0, nrow = n_cells, ncol = length(required_genes))
      colnames(X) <- required_genes
      rownames(X) <- colnames(binary_matrix)
      
      # Fill in available genes
      if (length(present_genes) > 0) {
        X[, present_genes] <- t(binary_matrix[present_genes, , drop = FALSE])
      }
    }
  } else if (method == "random_forest") {
    # For Random Forest: use sparse matrices for efficiency, convert to data.frame at end
    # (ranger requires data.frame, but we can keep things sparse until then)
    if (requireNamespace("Matrix", quietly = TRUE)) {
      # Create sparse matrix (more memory efficient for binary data)
      X <- Matrix::Matrix(0, nrow = n_cells, ncol = length(required_genes), sparse = TRUE)
      colnames(X) <- required_genes
      rownames(X) <- colnames(binary_matrix)
      
      # Fill in available genes using sparse matrix operations
      if (length(present_genes) > 0) {
        if (inherits(binary_matrix, "sparseMatrix")) {
          X[, present_genes] <- Matrix::t(binary_matrix[present_genes, , drop = FALSE])
        } else {
          X[, present_genes] <- Matrix::t(Matrix::Matrix(binary_matrix[present_genes, , drop = FALSE], sparse = TRUE))
        }
      }
      
      # Convert sparse matrix to dense data.frame for ranger (only at the end)
      X <- as.data.frame(as.matrix(X))
    } else {
      # Fallback to dense matrix if Matrix package not available
      X <- matrix(0, nrow = n_cells, ncol = length(required_genes))
      colnames(X) <- required_genes
      rownames(X) <- colnames(binary_matrix)
      
      # Fill in available genes
      if (length(present_genes) > 0) {
        X[, present_genes] <- t(as.matrix(binary_matrix[present_genes, , drop = FALSE]))
      }
      
      # Convert to data.frame for ranger
      X <- as.data.frame(X)
    }
  } else if (method == "xgboost") {
    # For XGBoost: use dense matrix (xgboost works better with dense matrices)
    X <- matrix(0, nrow = n_cells, ncol = length(required_genes))
    colnames(X) <- required_genes
    rownames(X) <- colnames(binary_matrix)
    
    # Fill in available genes
    if (length(present_genes) > 0) {
      # Convert to dense if sparse
      if (requireNamespace("Matrix", quietly = TRUE) && inherits(binary_matrix, "sparseMatrix")) {
        X[, present_genes] <- t(as.matrix(binary_matrix[present_genes, , drop = FALSE]))
      } else {
        X[, present_genes] <- t(binary_matrix[present_genes, , drop = FALSE])
      }
    }
  }
  # Missing genes remain 0 (default for binary: not expressed)
  
  if (verbose) {
    cat("Predicting position of", n_cells, "cells using ML model...\n")
    cat("  Model method:", method, "\n")
    cat("  Genes used:", trained_model$genes_used, "\n")
    cat("  Features:", length(required_genes), "genes\n")
  }
  
  # Make predictions using the trained model (method-specific)
  if (method == "ridge_regression") {
    # Check if glmnet is available
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("glmnet package is required for Ridge regression. Please install it using: install.packages('glmnet')")
    }
    
    predicted_pseudotime <- as.numeric(
      predict(trained_model$model,
              newx = X,
              s = "lambda.min")
    )
    
  } else if (method == "random_forest") {
    # Check if ranger is available
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("ranger package is required for Random Forest. Please install it using: install.packages('ranger')")
    }
    
    # Clean column names to match what was used during training
    # (make.names is deterministic, so same names will be cleaned the same way)
    colnames(X) <- make.names(colnames(X), unique = TRUE)
    
    # ranger prediction
    pred_result <- predict(trained_model$model, data = X)
    predicted_pseudotime <- as.numeric(pred_result$predictions)
    
  } else if (method == "xgboost") {
    # Check if xgboost is available
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("xgboost package is required for XGBoost. Please install it using: install.packages('xgboost')")
    }
    
    # xgboost prediction (requires DMatrix)
    dtest <- xgboost::xgb.DMatrix(data = X)
    predicted_pseudotime <- as.numeric(
      predict(trained_model$model, newdata = dtest)
    )
  } else {
    stop("Unknown method: ", method)
  }
  
  # Get cell names
  cell_names <- rownames(X)
  
  # Convert predictions to PPR_OBJECT format
  ppr_obj <- create_ppr_from_predictions(
    predicted_pseudotime = predicted_pseudotime,
    cell_names = cell_names,
    pseudotime_range = trained_model$pseudotime_range
  )
  
  if (verbose) {
    cat("Prediction complete.\n")
  }
  
  return(ppr_obj)
}
