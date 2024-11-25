#' @title reduce_counts_matrix
#'
#' @description
#' Reduces the counts matrix to only include genes included in switching_genes.
#'
#' @param sample_sce A Single Cell Experiment object,
#' containing a matrix of your samples gene expression.
#' @param switching_genes Genes which switch through the trajectory,
#' as identified by GeneSwitches.
#'
#' @return sample_sce A Single Cell Experiment object,
#' containing a reduced matrix of your samples gene expression,
#' only genes included in switching_gene are included.
#' 
#' @export
#'
reduce_counts_matrix <- function(sample_sce, switching_genes) {
  switching_genes_idx <- which(rownames(sample_sce) %in% switching_genes[, 1])
  sample_sce <- sample_sce[switching_genes_idx, , drop = FALSE]
  return(sample_sce)
}
