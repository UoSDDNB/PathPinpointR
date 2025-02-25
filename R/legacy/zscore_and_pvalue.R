#' @title Legacy Function: zscore_and_pvalue
#'
#' @description
#' Calculate the zscore and the p-value of the raw ppr score,
#'  by comparing to 2000 random samples.
#'
#' @param sce A Single Cell Experiment object,
#'  containing a matrix of your samples binary gene expression,
#'  which has been filtered to only include switching genes,
#'  using PathPinpointR::subset_switching_genes().
#' @param ppr An object of class ppr, must be the same sample.
#' @param switching_genes Genes which switch through the trajectory,
#' @param cpu Number of cores to use for parallel processing,
#'  (default is one less than available cores, 4 works best).
#'
#' @return the zscore_and_pvalue of the samples raw ppr score.
#'
#' @importFrom parallel clusterExport makeCluster parLapply stopCluster detectCores
#' @importFrom stats sd
#' @importFrom methods is
#' @export
#'

zscore_and_pvalue <- function(sce, ppr, switching_genes, cpu = 1) {
  usethis::deprecated("This is a legacy function and will not be maintained actively.")
  
  ## checks

  # check that sce is a SingleCellExperiment object
  if (!is(sce, "SingleCellExperiment")) {
    stop("\n  The sce argument must be a SingleCellExperiment object.")
  }

  # Check if sce contains binarized expression data
  if (is.null(sce@assays@data@listData$binary)) {
    # If binarized expression data is missing, display a message and stop.
    stop("\n  Binarized expression data not found in \"",
         deparse(substitute(sce)),
         "\" \n  You must run `GeneSwitches::binarize_exp()` first.")
  }

  # check that sce has been subset
  if (dim(sce)[1] > dim(switching_genes)[1]) {
    stop("\n  The number of genes in the reduced binary counts matrix is
     greater than the number of switching genes.
     
    Make sure you have run PathPinpointR::subset_switching_genes()")
  }

  # check that ppr is a ppr object
  if (!is(ppr, "PPR_OBJECT")) {
    stop("\n  The ppr argument must be a PPR_OBJECT.")
  }

  # check that sce and ppr come from the same sample
  if (dim(sce)[1] != dim(ppr[[1]][[1]])[1] || dim(sce)[2] != length(ppr[[1]])) {
    stop("\n  The sce and ppr objects must be from the same sample.")
  }

  # check that switching_genes is a > class(switching_genes)
  if (!is(switching_genes, "DFrame")) {
    cat("Warning:
      Make sure switching_genes argument is a DFrame.\n")
  }

  # if cpu is set to 0 or is larger or equal to than the number available.
  # use all but one available cores
  if (cpu == 0 || cpu >= detectCores()) {
    cpu <- detectCores() - 1
  }

  # Find the maximum raw ppr score.
  max_raw_ppr_score <- max(ppr$sample_flat)

  # Extract the binary matrix
  bin_mat <- sce@assays@data$binary

  # Define the number of iterations
  n_iterations <- 2000

  # Create a cluster for parallel processing
  cl <- makeCluster(cpu)

  # Cluster setup
  clusterExport(cl,
                c("predict_position",
                  "switching_genes",
                  "which_mid_max",
                  "sce"),
                  envir = environment())

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

  # assign Z-score to the ppr object
  ppr$z_score <- z_score
  
  # calculate how many values in random_max_raw_ppr_scores,
  # are greater than max_raw_ppr_score
  p_value <- sum(random_max_raw_ppr_scores > max_raw_ppr_score) / n_iterations
  # if p_value is 0, set it to "<0.0005"
  if (p_value == 0) {
    p_value <- "<0.0005"
  }
  # assign p-value to the ppr object
  ppr$p_value <- p_value

  #return the ppr object
  return(ppr)
}