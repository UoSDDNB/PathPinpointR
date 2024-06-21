#' @title zscore
#'
#' @description
#' Calculate the zscore of the raw ppr score,
#'  by comparing to 2000 random samples.
#'
#' @param sce A Single Cell Experiment object,
#'  containing a matrix of your samples binary gene expression,
#'  which has been filtered to only include switching genes,
#'  using PathPinPointR::subset_switching_genes().
#' @param ppr An object of class ppr, must be the same sample.
#' @param switching_genes Genes which switch through the trajectory,
#' @param cpu Number of cores to use for parallel processing,
#'  (default is one less than available cores).
#'
#' @return the zscore of the samples raw ppr score.
#'
#' @importFrom parallel makeCluster parLapply stopCluster
#' @export
#'

calculate_zscore <- function(sce, ppr, switching_genes, cpu = 0) {

  ## check

  # check that sce is a SingleCellExperiment object

  # check that sce is reduced?

  # check that ppr is a ppr object

  # check that sce and ppr come from the same sample

  # check that switching_genes is a df?

  # if cpu is set to 0 or is larger or equal to than the number available.
  # use all but one available cores
  if (cpu == 0 || cpu >= detectCores()) {
    cpu <- paralell::detectCores() - 1
  }


  # Find the maximum raw ppr score.
  max_raw_ppr_score <- max(ppr$sample_flat)

  # Extract the binary matrix
  bin_mat <- sce@assays@data$binary

  

  # Define the number of iterations
  n_iterations <- 2000

  # Create a cluster for parallel processing
  cl <- parallel::makeCluster(cpu)

  # Cluster setup
  parallel::clusterExport(cl,
                          c("predict_position",
                            "switching_genes",
                            "which_mid_max",
                            "sce"))

  # Parallel loop using parLapply
  random_max_raw_ppr_scores <- parLapply(cl, 1:n_iterations, function(i) {

    # Shuffle the row names in place
    bin_mat <- sce@assays@data$binary
    rownames(bin_mat) <- sample(rownames(bin_mat))
    sce@assays@data$binary <- bin_mat
    # Perform prediction
    random_ppr <- predict_position(sce, switching_genes)
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