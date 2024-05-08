#' Subset Gene Expression for Switching Genes
#'
#' @description
#' Reduce a binary expression matrix to only the switching genes.
#'
#' @param sample_sce a Single Cell Experiment object,
#' which contains a binarized expression matrix,
#' this is your sample.
#'
#' @param switching_genes Genes identified by GeneSwitches as "switching",
#' produced using your reference.
#'
#' @return a binary expression matrix subseted to only include switching genes
#' @export
#'
#'

subset_switching_genes <- function(sample_sce, switching_genes) {
  # Check if sample_sce contains binarized expression data
  if (is.null(sample_sce@assays@data@listData$binary)) {
    # If binarized expression data is missing, display a message and stop.
    stop("\n  Binarized expression data not found in \"",
         deparse(substitute(sample_sce)),
         "\" \n  You must run `GeneSwitches::binarize_exp()` first.")
  }
  switching_genes_idx <- which(rownames(sample_sce) %in% switching_genes[, 1])
  reduced_sample_sce <- sample_sce[switching_genes_idx, , drop = FALSE]
  return(reduced_sample_sce)
}
