# Test script for ML-based prediction approach
# This script follows the PathPinpointR tutorial exactly and adds ML approach testing

# Load required packages
library(PathPinpointR)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(GeneSwitches)
library(ggplot2)
library(RColorBrewer)

## Load the reference data

# download the example data
get_example_data()

# Load the reference data to the environment
reference_seu <- readRDS("./reference.rds")
# Load the sample data to the environment
samples_seu <- lapply(c("./LW120.rds", "./LW122.rds"), readRDS)
# name the samples in the list
names(samples_seu) <- c("LW120", "LW122")

#### View the reference UMAP plot

# Plot the reference data, colored by the day of development.
DimPlot(object = reference_seu,
        reduction = "umap",
        group.by = "time",
        label = TRUE) +
  ggtitle("Reference")

## Convert to SingleCellExperiment objects.

# Prior to running slingshot and GeneSwitches, we need to convert the
# Seurat objects to SingleCellExperiment objects.
reference_sce <- SingleCellExperiment(assays = list(expdata = reference_seu@assays$RNA$counts))
colData(reference_sce) <- DataFrame(reference_seu@meta.data)
reducedDims(reference_sce)$UMAP <- reference_seu@reductions$umap@cell.embeddings

# create an empty list to to store the sce sample objects
samples_sce <- list()
# Iterate through each Seurat object in the samples list
for (i in seq_along(samples_seu)){
  # convert each sample to a SingleCellExperiment object & store in the list
  samples_sce[[i]] <- SingleCellExperiment(assays = list(expdata = samples_seu[[i]]@assays$RNA$counts))
}
# carry over the sample names from the Seurat objects.
names(samples_sce) <- names(samples_seu)

## Run slingshot

# Run slingshot on the reference data to produce pseudotime for each cell.
reference_sce <- slingshot(reference_sce,
                            clusterLabels = "seurat_clusters",
                            start.clus  = "2",
                            end.clus = "1",
                            reducedDim = "UMAP")

#Rename the Pseudotime column to work with GeneSwitches
colData(reference_sce)$Pseudotime <- reference_sce$slingPseudotime_1

#### Plot the slingshot trajectory.

# The plot shows the trajectory of the blastocyst data, with cells colored
# by pseudotime.
# Generate colors
colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)
plotcol <- colors[cut(reference_sce$slingPseudotime_1, breaks = 100)]
# Plot the data
plot(reducedDims(reference_sce)$UMAP, col = plotcol, pch = 16, asp = 1)
lines(SlingshotDataSet(reference_sce), lwd = 2, col = "black")

## Binarize the Gene Expression Data

# Using the package GeneSwitches, binarize the gene expression data of the 
# reference and query data-sets, with a cutoff of 1.
# binarize the expression data of the reference
reference_sce <- binarize_exp(reference_sce,
                              fix_cutoff = TRUE,
                              binarize_cutoff = 1,
                              ncores = 1)

# Find the switching point of each gene in the reference data
reference_sce <- find_switch_logistic_fastglm(reference_sce,
                                              downsample = TRUE)

## Visualise the switching genes

# generate a list of switching genes, and visualise them on a pseudo timeline.
switching_genes <- filter_switchgenes(reference_sce,
                                      allgenes = TRUE,
                                      r2cutoff = 0)

# Check if switching genes were found
cat("Number of switching genes found:", nrow(switching_genes), "\n")
if (nrow(switching_genes) == 0) {
  warning("No switching genes found! This may cause issues later.")
  cat("Trying with more permissive settings...\n")
  switching_genes <- filter_switchgenes(reference_sce,
                                        allgenes = TRUE,
                                        r2cutoff = -10,
                                        topnum = NULL)
  cat("Number of switching genes with permissive settings:", nrow(switching_genes), "\n")
}

# Plot the timeline using plot_timeline_ggplot
plot_timeline_ggplot(switching_genes,
                     timedata = colData(reference_sce)$Pseudotime,
                     txtsize = 3)

# Individual switching genes can be visualised using Seurat's FeaturePlot function.
FeaturePlot(reference_seu, 
            features = rownames(switching_genes)[1], 
            reduction = "umap", 
            min.cutoff = "q10", 
            max.cutoff = "q90") + 
  ggtitle(rownames(switching_genes)[1])

## Select a number of switching genes

# Using the PPR function precision():
precision(reference_sce)

# Narrow down the search to find the optimum number of switching genes.
precision(reference_sce, n_sg_range = seq(50, 150, 1))

## Filter and re-visualise the switching genes.

# Filter to top switching genes (using example value from tutorial)
# Only filter if we have switching genes
if (nrow(switching_genes) > 0) {
  switching_genes <- filter_switchgenes(reference_sce,
                                        allgenes = TRUE,
                                        r2cutoff = 0,
                                        topnum = 114)
  cat("Number of switching genes after filtering to top 114:", nrow(switching_genes), "\n")
} else {
  warning("switching_genes is empty - cannot filter to top 114")
}

# Plot the timeline using plot_timeline_ggplot
plot_timeline_ggplot(switching_genes,
                     timedata = colData(reference_sce)$Pseudotime,
                     txtsize = 3)

## Measure accuracy (using current method)

# reduce the reference data to only include the switching genes.
reference_reduced_sce <- reduce_counts_matrix(reference_sce, switching_genes)
# predict the position of cells in the reference trajectory.
reference_ppr <- predict_position(reference_reduced_sce, switching_genes)

# plot the accuracy of the prediction
accuracy_test(reference_ppr, reference_sce, plot = TRUE)

# ============================================================================
# ML APPROACH TESTING
# ============================================================================

cat("\n=== ML APPROACH TESTING ===\n\n")

## Split reference data into train/test sets (80/20)
set.seed(123)
train_idx <- sample(seq_len(ncol(reference_sce)), size = 0.8 * ncol(reference_sce))
train_sce <- reference_sce[, train_idx]
test_sce <- reference_sce[, -train_idx]

cat("Training set:", length(train_idx), "cells\n")
cat("Test set:", ncol(test_sce), "cells\n\n")

## Debug: Check switching_genes before training
cat("=== DEBUGGING SWITCHING_GENES ===\n")
cat("switching_genes is null:", is.null(switching_genes), "\n")
if (!is.null(switching_genes)) {
  cat("switching_genes class:", paste(class(switching_genes), collapse = ", "), "\n")
  cat("switching_genes is data.frame:", is.data.frame(switching_genes), "\n")
  cat("switching_genes is DFrame:", inherits(switching_genes, "DFrame"), "\n")
  # Check if it's a data.frame-like object
  is_df_like <- is.data.frame(switching_genes) || 
                inherits(switching_genes, "DFrame") || 
                inherits(switching_genes, "DataFrame")
  if (is_df_like) {
    cat("switching_genes nrow:", nrow(switching_genes), "\n")
    cat("switching_genes ncol:", ncol(switching_genes), "\n")
    if (nrow(switching_genes) > 0) {
      cat("First few rownames:", head(rownames(switching_genes)), "\n")
    }
  }
}
cat("===============================\n\n")

## Train ML models on training set (all three methods)
cat("Training ML models (Ridge, Random Forest, XGBoost)...\n\n")

# Check if switching_genes is valid, otherwise use all genes
# DFrame objects are not data.frames but work similarly
is_df_like <- !is.null(switching_genes) && 
              (is.data.frame(switching_genes) || 
               inherits(switching_genes, "DFrame") || 
               inherits(switching_genes, "DataFrame"))

# Prepare training arguments
if (!is_df_like || nrow(switching_genes) == 0) {
  cat("WARNING: switching_genes is empty or invalid. Using all genes instead.\n")
  train_args <- list(genes = "all")
} else {
  train_args <- list(switching_genes = switching_genes, genes = "switching")
}

# Train Ridge regression model
cat("=== Training Ridge Regression ===\n")
trained_model_ridge <- do.call(train_ml_model, 
                               c(list(reference_sce = train_sce,
                                      method = "ridge",
                                      alpha = 0,
                                      cv_folds = 5),
                                 train_args))
print(trained_model_ridge)
cat("\n")

# Train Random Forest model
cat("=== Training Random Forest ===\n")
trained_model_rf <- do.call(train_ml_model,
                            c(list(reference_sce = train_sce,
                                   method = "random_forest",
                                   num.trees = 500,
                                   mtry = NULL,  # Will use default sqrt(n_features)
                                   min.node.size = 5),
                              train_args))
print(trained_model_rf)
cat("\n")

# Train XGBoost model
cat("=== Training XGBoost ===\n")
trained_model_xgb <- do.call(train_ml_model,
                             c(list(reference_sce = train_sce,
                                    method = "xgboost",
                                    nrounds = 100,
                                    max_depth = 6,
                                    eta = 0.3,
                                    subsample = 0.8),
                               train_args))
print(trained_model_xgb)
cat("\n")

## Predict on test set using all ML methods
cat("=== PREDICTING WITH ML MODELS ===\n\n")

cat("Predicting with Ridge regression...\n")
test_ppr_ridge <- predict_position_ml(test_sce, trained_model_ridge)

cat("Predicting with Random Forest...\n")
test_ppr_rf <- predict_position_ml(test_sce, trained_model_rf)

cat("Predicting with XGBoost...\n")
test_ppr_xgb <- predict_position_ml(test_sce, trained_model_xgb)

## Predict on test set using current method (for comparison)
cat("\nPredicting with current method...\n")
test_reduced_sce <- reduce_counts_matrix(test_sce, switching_genes)
test_ppr_current <- predict_position(test_reduced_sce, switching_genes)

## Compare accuracy
cat("\n=== ACCURACY COMPARISON ===\n\n")

# Ridge regression accuracy
cat("Ridge Regression:\n")
ridge_accuracy <- accuracy_test(test_ppr_ridge, test_sce, plot = TRUE)

# Random Forest accuracy
cat("\nRandom Forest:\n")
rf_accuracy <- accuracy_test(test_ppr_rf, test_sce, plot = TRUE)

# XGBoost accuracy
cat("\nXGBoost:\n")
xgb_accuracy <- accuracy_test(test_ppr_xgb, test_sce, plot = TRUE)

# Current method accuracy
cat("\nCurrent Method:\n")
current_accuracy <- accuracy_test(test_ppr_current, test_sce, plot = TRUE)

# Simple comparison message
cat("\n=== COMPARISON SUMMARY ===\n")
cat("All ML methods have been evaluated. Compare the plots above to see the differences.\n")
cat("Lower inaccuracy values (shown in the plots) indicate better predictions.\n")
cat("\nMethods tested:\n")
cat("  1. Ridge Regression (default ML method)\n")
cat("  2. Random Forest\n")
cat("  3. XGBoost\n")
cat("  4. Current method (original PathPinpointR approach)\n")

## Test on independent sample data

cat("\n=== TESTING ON INDEPENDENT SAMPLES ===\n\n")

# Binarize the sample data
# First reduce the sample data to only include the switching genes.
samples_sce <- lapply(samples_sce, reduce_counts_matrix, switching_genes)

# binarize the expression data of the samples
samples_binarized <- lapply(samples_sce,
                            binarize_exp,
                            fix_cutoff = TRUE,
                            binarize_cutoff = 1,
                            ncores = 1)

# Predict using ML approaches (using Ridge as default, but can test others)
cat("Predicting sample positions with ML models...\n")
cat("Using Ridge regression for sample predictions (can change to trained_model_rf or trained_model_xgb)\n")
samples_ppr_ridge <- lapply(samples_binarized, 
                           predict_position_ml, 
                           trained_model = trained_model_ridge)
samples_ppr_rf <- lapply(samples_binarized, 
                        predict_position_ml, 
                        trained_model = trained_model_rf)
samples_ppr_xgb <- lapply(samples_binarized, 
                         predict_position_ml, 
                         trained_model = trained_model_xgb)

# Predict using current method
cat("Predicting sample positions with current method...\n")
samples_ppr_current <- lapply(samples_binarized, 
                              predict_position, 
                              switching_genes = switching_genes)

## Plotting the predicted position of each sample

# Ridge regression predictions
cat("\nPlotting Ridge regression predictions...\n")
ppr_plot() +
  reference_idents(reference_sce, "time") +
  sample_prediction(samples_ppr_ridge[[1]], label = "Ridge: LW120", col = "red") +
  sample_prediction(samples_ppr_ridge[[2]], label = "Ridge: LW122", col = "orange")

# Random Forest predictions
cat("Plotting Random Forest predictions...\n")
ppr_plot() +
  reference_idents(reference_sce, "time") +
  sample_prediction(samples_ppr_rf[[1]], label = "RF: LW120", col = "purple") +
  sample_prediction(samples_ppr_rf[[2]], label = "RF: LW122", col = "magenta")

# XGBoost predictions
cat("Plotting XGBoost predictions...\n")
ppr_plot() +
  reference_idents(reference_sce, "time") +
  sample_prediction(samples_ppr_xgb[[1]], label = "XGB: LW120", col = "darkgreen") +
  sample_prediction(samples_ppr_xgb[[2]], label = "XGB: LW122", col = "green")

# Current method predictions
cat("Plotting current method predictions...\n")
ppr_plot() +
  reference_idents(reference_sce, "time") +
  sample_prediction(samples_ppr_current[[1]], label = "Current: LW120", col = "blue") +
  sample_prediction(samples_ppr_current[[2]], label = "Current: LW122", col = "cyan")

# Combined comparison plot (Ridge vs Current)
cat("Plotting combined comparison (Ridge vs Current)...\n")
ppr_plot() +
  reference_idents(reference_sce, "time") +
  sample_prediction(samples_ppr_ridge[[1]], label = "Ridge: LW120", col = "red") +
  sample_prediction(samples_ppr_current[[1]], label = "Current: LW120", col = "blue") +
  switching_times(c(rownames(switching_genes)[1:3]), switching_genes)

# All ML methods comparison
cat("Plotting all ML methods comparison...\n")
ppr_plot() +
  reference_idents(reference_sce, "time") +
  sample_prediction(samples_ppr_ridge[[1]], label = "Ridge: LW120", col = "red") +
  sample_prediction(samples_ppr_rf[[1]], label = "RF: LW120", col = "purple") +
  sample_prediction(samples_ppr_xgb[[1]], label = "XGB: LW120", col = "darkgreen")

# Vioplot comparisons
cat("Plotting vioplot comparisons...\n")
ppr_vioplot(samples_ppr_ridge, reference_sce, ident = "time")
title("Ridge Regression Predictions")

ppr_vioplot(samples_ppr_rf, reference_sce, ident = "time")
title("Random Forest Predictions")

ppr_vioplot(samples_ppr_xgb, reference_sce, ident = "time")
title("XGBoost Predictions")

ppr_vioplot(samples_ppr_current, reference_sce, ident = "time")
title("Current Method Predictions")

cat("\n=== TESTING COMPLETE ===\n")



