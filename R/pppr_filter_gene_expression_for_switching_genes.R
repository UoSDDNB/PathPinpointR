#' Filter Gene Expression for Switching Genes
#'
#' @description
#' Create a reduced binary expression matrix for only the selected switching genes,
#'     binary_counts_matrix is from the sample DATA and gs_scorer_genes is from Atlas Data.
#'
#' @param binary_counts_matrix a binary expression matrix from your sample.
#' @param gs_scorer_genes Switching genes which are evenly distributed through pseudotime..
#'
#' @return a reduced binary expression matrix filtered to only include selected switching genes
#' @export
#'
#'

filter_gene_expression_for_switching_genes <- function(binary_counts_matrix, reference.sg) {
  indices_of_switching_genes   <- which(rownames(binary_counts_matrix) %in% reference.sg[,1])
  reduced_binary_counts_matrix <- binary_counts_matrix[indices_of_switching_genes, ,drop = FALSE]
  return(reduced_binary_counts_matrix)
}
