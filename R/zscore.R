#' @title zscore
#'
#' @description
#' Calculate the zscore of the raw ppr score.
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

calculate_zscore <- function(sce, ppr, switching_genes) {
  # Find the maximum raw ppr score.
  max_raw_ppr_score <- max(ppr$sample_flat)

  ## produce 2000 random samples, and extract their max ppr scores.
  # Extract the binary matrix
  bin_mat <- sce@assays@data$binary

  # Create a vector to store the max ppr scores of the random samples
  random_max_raw_ppr_scores <- numeric(2000)

  # predict the position of 2000 random samples
  # and extract the max ppr score
  for (i in 1:2000) {
    # Shuffle the row names in place
    rownames(bin_mat) <- sample(rownames(bin_mat))

    # Update the binary matrix within the sample object
    sce@assays@data$binary <- bin_mat

    # Perform prediction
    random_ppr <- ppr_predict_position(sce, switching_genes)

    # Extract the maximum raw ppr score
    random_max_raw_ppr_scores[i] <- max(random_ppr$sample_flat)
  }

  # calculate the distance of the max_raw_ppr_score from the mean of random_max_raw_ppr_scores
  # using standard deviations of the random samples

  # calculate the standand deviation of the random max ppr scores
  sd_random_max_raw_ppr_scores <- sd(random_max_raw_ppr_scores)

  # calculate the z-score
  z_score <- (max_raw_ppr_score - mean(random_max_raw_ppr_scores)) / sd_random_max_raw_ppr_scores

  return(z_score)
}