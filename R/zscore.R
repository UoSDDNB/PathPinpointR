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
#' @param cpu The number of cores to use for parallel processing.
#'
#' @return the zscore of the raw ppr score.
#'
#' @importFrom stats shapiro.test
#' @importFrom parallel makeCluster parLapply stopCluster
#' @export
#'

calculate_zscore <- function(sce, ppr, switching_genes, cpu = 0) {
  # Find the maximum raw ppr score.
  max_raw_ppr_score <- max(ppr$sample_flat)

  ## produce 2000 random samples, and extract their max ppr scores.
  # Extract the binary matrix
  bin_mat <- sce@assays@data$binary

  # if cpu is set to 0 use all but one available cores
  if (cpu == 0) {
    cpu <- detectCores() - 1
  }

  # Define the number of iterations
  n_iterations <- 2000

  # Create a cluster for parallel processing
  cl <- makeCluster(cpu)

  # Cluster setup
  clusterExport(cl,
                c("predict_position",
                  "switching_genes",
                  "sce",
                  "which_mid_max"))

  # Parallel loop using parLapply
  random_max_raw_ppr_scores <- parLapply(cl, 1:n_iterations, function(i) {

    # Copy the SingleCellExperiment object for each iteration
    sce_copy <- sce
    # Shuffle the row names in place
    bin_mat <- sce_copy@assays@data$binary
    rownames(bin_mat) <- sample(rownames(bin_mat))
    sce_copy@assays@data$binary <- bin_mat
    # Perform prediction
    random_ppr <- predict_position(sce_copy, switching_genes)
    # Extract the maximum raw ppr score
    max(random_ppr$sample_flat)
  })

  # Stop the cluster
  stopCluster(cl)

  # unlist the random_max_raw_ppr_scores (make it a vector)
  random_max_raw_ppr_scores <- unlist(random_max_raw_ppr_scores)

  # calculate the distance of the max_raw_ppr_score from the mean of random_max_raw_ppr_scores
  # using standard deviations of the random samples

  # calculate the standard deviation of the random max ppr scores
  sd_random_max_raw_ppr_scores <- sd(random_max_raw_ppr_scores)

  # calculate the z-score
  z_score <- (max_raw_ppr_score - mean(random_max_raw_ppr_scores)) / sd_random_max_raw_ppr_scores

  return(z_score)
}