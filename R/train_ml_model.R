#' @title Train ML Model for Pseudotime Prediction
#'
#' @description
#' Trains a machine learning model to predict pseudotime from binary gene expression.
#' Supports three methods: Ridge regression, Random Forest, and XGBoost.
#' The model can use either switching genes only or all available genes.
#'
#' @param reference_sce A SingleCellExperiment object containing:
#' \itemize{
#'   \item Binarized expression data (from \code{GeneSwitches::binarize_exp()})
#'   \item Pseudotime information (from slingshot or similar)
#' }
#' @param switching_genes Optional. Data frame of switching genes (as output from
#' \code{GeneSwitches::filter_switchgenes()}). Required if \code{genes = "switching"}.
#' @param genes Character. Which genes to use as features:
#' \itemize{
#'   \item \code{"switching"} (default): Use only switching genes
#'   \item \code{"all"}: Use all genes in the binary expression matrix
#' }
#' @param method Character. ML method to use:
#' \itemize{
#'   \item \code{"ridge"} (default): Ridge regression using glmnet
#'   \item \code{"random_forest"}: Random Forest using ranger
#'   \item \code{"xgboost"}: Gradient boosting using xgboost
#' }
#' @param alpha Numeric. Elastic net mixing parameter (for Ridge regression only):
#' \itemize{
#'   \item \code{0} (default): Ridge regression
#'   \item \code{1}: Lasso regression
#'   \item \code{0 < alpha < 1}: Elastic net
#' }
#' @param cv_folds Integer. Number of folds for cross-validation (default: 5, for Ridge only).
#' @param ... Additional hyperparameters for method-specific training:
#' \itemize{
#'   \item For Random Forest: \code{num.trees} (default: 500), \code{mtry} (default: sqrt(n_features)), \code{min.node.size} (default: 5)
#'   \item For XGBoost: \code{nrounds} (default: 100), \code{max_depth} (default: 6), \code{eta} (default: 0.3), \code{subsample} (default: 0.8)
#' }
#'
#' @return A trained model object (list) containing:
#' \itemize{
#'   \item \code{model}: The trained model object (varies by method)
#'   \item \code{gene_features}: Character vector of gene names used as features
#'   \item \code{pseudotime_range}: Numeric vector c(min, max) of training pseudotime range
#'   \item \code{pseudotime_step}: Calculated step size for index conversion
#'   \item \code{method}: Character string ("ridge_regression", "random_forest", or "xgboost")
#'   \item \code{genes_used}: Character string indicating which genes were used
#' }
#'
#' @details
#' This function trains a machine learning model to predict continuous pseudotime values
#' from binary gene expression patterns. Three methods are supported:
#' \itemize{
#'   \item \strong{Ridge regression}: Uses cross-validation to select optimal regularization parameter
#'   \item \strong{Random Forest}: Ensemble of decision trees for robust predictions
#'   \item \strong{XGBoost}: Gradient boosting for high-performance predictions
#' }
#'
#' @examples
#' \dontrun{
#' # Train Ridge regression (default) on reference data with switching genes
#' trained_model <- train_ml_model(reference_sce, 
#'                                  switching_genes = switching_genes,
#'                                  genes = "switching")
#'
#' # Train Random Forest
#' model_rf <- train_ml_model(reference_sce, 
#'                            switching_genes = switching_genes,
#'                            method = "random_forest")
#'
#' # Train XGBoost with custom hyperparameters
#' model_xgb <- train_ml_model(reference_sce,
#'                             switching_genes = switching_genes,
#'                             method = "xgboost",
#'                             nrounds = 200,
#'                             max_depth = 8)
#' }
#'
#' @importFrom glmnet cv.glmnet
#' @export
train_ml_model <- function(reference_sce,
                           switching_genes = NULL,
                           genes = c("switching", "all"),
                           method = c("ridge", "random_forest", "xgboost"),
                           alpha = 0,
                           cv_folds = 5,
                           ...) {
  
  # Validate method parameter
  method <- match.arg(method)
  
  # Validate genes parameter
  genes <- match.arg(genes)
  
  # Check for switching genes if needed
  if (genes == "switching") {
    if (is.null(switching_genes)) {
      stop("switching_genes must be provided when genes = 'switching'")
    }
    # Extract gene names from switching_genes
    # switching_genes can be a data.frame or DFrame (Bioconductor) with gene names as rownames
    # Check if it's a data.frame-like object (data.frame, DFrame, DataFrame, etc.)
    is_df_like <- is.data.frame(switching_genes) || 
                  inherits(switching_genes, "DFrame") || 
                  inherits(switching_genes, "DataFrame")
    
    if (is_df_like && nrow(switching_genes) > 0) {
      gene_features <- rownames(switching_genes)
    } else {
      stop("switching_genes must be a non-empty data frame (or DFrame) with gene names as rownames")
    }
  } else {
    gene_features <- NULL  # Will use all genes
  }
  
  # Prepare training data
  training_data <- prepare_training_data(reference_sce, gene_features)
  X <- training_data$X
  y <- training_data$y
  
  # Store pseudotime range for later use
  pseudotime_range <- c(min(y), max(y))
  pseudotime_step <- (pseudotime_range[2] - pseudotime_range[1]) / 100
  
  # Get final gene features used
  final_gene_features <- colnames(X)
  
  # Get additional hyperparameters from ...
  hyperparams <- list(...)
  
  # Train model based on method
  if (method == "ridge") {
    # Check if glmnet is available
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("glmnet package is required for Ridge regression. Please install it using: install.packages('glmnet')")
    }
    
    # Train Ridge regression model with cross-validation
    cat("Training Ridge regression model...\n")
    cat("  Features:", length(final_gene_features), "genes\n")
    cat("  Samples:", nrow(X), "cells\n")
    cat("  Cross-validation folds:", cv_folds, "\n")
    
    # Train model
    model <- glmnet::cv.glmnet(x = X,
                               y = y,
                               alpha = alpha,
                               nfolds = cv_folds,
                               standardize = TRUE)
    
    cat("  Optimal lambda:", model$lambda.min, "\n")
    cat("  Model training complete.\n")
    
    method_name <- "ridge_regression"
    
  } else if (method == "random_forest") {
    # Check if ranger is available
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("ranger package is required for Random Forest. Please install it using: install.packages('ranger')")
    }
    
    # Convert to data.frame for ranger (ranger requires data.frame input)
    # Convert sparse matrix to dense if needed
    if (requireNamespace("Matrix", quietly = TRUE) && inherits(X, "sparseMatrix")) {
      X_df <- as.data.frame(as.matrix(X))
    } else {
      X_df <- as.data.frame(X)
    }
    
    # Set default hyperparameters
    num_trees <- if (is.null(hyperparams$num.trees)) 500 else hyperparams$num.trees
    mtry <- if (is.null(hyperparams$mtry)) floor(sqrt(ncol(X_df))) else hyperparams$mtry
    min_node_size <- if (is.null(hyperparams$min.node.size)) 5 else hyperparams$min.node.size
    
    cat("Training Random Forest model...\n")
    cat("  Features:", length(final_gene_features), "genes\n")
    cat("  Samples:", nrow(X_df), "cells\n")
    cat("  Number of trees:", num_trees, "\n")
    cat("  mtry:", mtry, "\n")
    cat("  Min node size:", min_node_size, "\n")
    
    # Clean column names to be valid R identifiers (avoids formula issues)
    colnames(X_df) <- make.names(colnames(X_df), unique = TRUE)
    
    # Create data frame with dependent variable
    train_data <- cbind(y = y, X_df)
    
    # Train model using ranger's alternative interface (no formula)
    # This avoids issues with special characters in gene names
    model <- ranger::ranger(
      dependent.variable.name = "y",
      data = train_data,
      num.trees = num_trees,
      mtry = mtry,
      min.node.size = min_node_size,
      verbose = FALSE
    )
    
    cat("  Model training complete.\n")
    cat("  OOB prediction error:", round(model$prediction.error, 4), "\n")
    
    method_name <- "random_forest"
    
  } else if (method == "xgboost") {
    # Check if xgboost is available
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("xgboost package is required for XGBoost. Please install it using: install.packages('xgboost')")
    }
    
    # Convert to dense matrix for xgboost (xgboost works better with dense matrices)
    if (requireNamespace("Matrix", quietly = TRUE) && inherits(X, "sparseMatrix")) {
      X_dense <- as.matrix(X)
    } else {
      X_dense <- as.matrix(X)
    }
    
    # Set default hyperparameters
    nrounds <- if (is.null(hyperparams$nrounds)) 100 else hyperparams$nrounds
    max_depth <- if (is.null(hyperparams$max_depth)) 6 else hyperparams$max_depth
    eta <- if (is.null(hyperparams$eta)) 0.3 else hyperparams$eta
    subsample <- if (is.null(hyperparams$subsample)) 0.8 else hyperparams$subsample
    
    cat("Training XGBoost model...\n")
    cat("  Features:", length(final_gene_features), "genes\n")
    cat("  Samples:", nrow(X_dense), "cells\n")
    cat("  Number of rounds:", nrounds, "\n")
    cat("  Max depth:", max_depth, "\n")
    cat("  Learning rate (eta):", eta, "\n")
    cat("  Subsample:", subsample, "\n")
    
    # Create DMatrix for xgboost
    dtrain <- xgboost::xgb.DMatrix(data = X_dense, label = y)
    
    # Set parameters
    params <- list(
      objective = "reg:squarederror",
      max_depth = max_depth,
      eta = eta,
      subsample = subsample,
      verbose = 0
    )
    
    # Train model
    model <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = nrounds,
      verbose = 0
    )
    
    cat("  Model training complete.\n")
    
    method_name <- "xgboost"
  }
  
  # Create model object
  trained_model <- list(
    model = model,
    gene_features = final_gene_features,
    pseudotime_range = pseudotime_range,
    pseudotime_step = pseudotime_step,
    method = method_name,
    genes_used = genes
  )
  
  class(trained_model) <- "PPR_ML_MODEL"
  
  return(trained_model)
}
