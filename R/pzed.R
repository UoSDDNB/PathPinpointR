#' @title pzed
#'
#' @description
#' Calculate the zscore of the raw ppr score. parallel
#'
#' @param sce A Single Cell Experiment object,
#' containing a matrix of your samples binary gene expression,
#' which has been filtered to only include switching genes,
#' using PathPinPointR::subset_switching_genes().
#' @param ppr An object of class ppr, must be the same sample.
#' @param switching_genes Genes which switch through the trajectory,
#'
#' @return the zscore of the raw ppr score.
#'
#' @export
#'
#'

# Load necessary packages
library(future)
library(future.apply)


# Set up a parallel backend using a specified number of workers (CPUs)
future::plan(future::multisession, workers = 3)

pzed <- function(sce, ppr, switching_genes) {
  # Find the maximum raw ppr score
  max_raw_ppr_score <- max(ppr$sample_flat)
  
  # Extract the binary matrix
  bin_mat <- sce@assays@data$binary
  
  # Define a function to process each random sample
  process_sample <- function(i) {
    # Shuffle the row names in place
    rownames(bin_mat) <- sample(rownames(bin_mat))
    
    # Update the binary matrix within the sample object
    sce@assays@data$binary <- bin_mat
    
    # Perform prediction
    random_ppr <- predict_position(sce, switching_genes)
    
    # Extract the maximum raw ppr score
    max(random_ppr$sample_flat)
  }
  
  # Use future_lapply to parallelize the sampling
  random_max_raw_ppr_scores <- future.apply::future_lapply(1:2000, process_sample)
  
  # Convert the result to numeric
  random_max_raw_ppr_scores <- unlist(random_max_raw_ppr_scores)
  
  # Calculate the standard deviation of the random max ppr scores
  sd_random_max_raw_ppr_scores <- sd(random_max_raw_ppr_scores)
  
  # Calculate the z-score
  z_score <- (max_raw_ppr_score - mean(random_max_raw_ppr_scores)) / sd_random_max_raw_ppr_scores
  
  return(z_score)
}
