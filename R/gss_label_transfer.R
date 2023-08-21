
#' Label Transfer
#' To make sure that your labels in your patient data match that used in the atlas.
#'
#' @import Seurat
#' @param query The Query object. e.g. the patient data
#' @param reference The reference Seurat object. e.g the atlas data
#'
#' @return The Query object with new predicted cell types.
#' @export
label_transfer <- function(query, reference) {
  ## Find transfer anchors
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:20
  )

  ## Map Query
  query <- MapQuery(
    anchorset = anchors,
    query = query,
    reference = reference,
    refdata = list(
      cell_type = "lv1_annot"
    ),
    reference.reduction = "pca"
  )

  ## Remove cells from the patient data which don't have a confident cell_type prediction.
  query <- subset(query, subset = query$predicted.cell_type.score > 0.65)

  ## Process the patient data using Seurat

  # I'm not sure if SCTransform would be appropriate here.
  # query <- SCTransform(query)

  # Run PCA on the patient data
  query <- RunPCA(query)

  # Find nearest neighbors in the patient data
  query <- FindNeighbors(query)

  # Cluster the patient data
  query <- FindClusters(query)

  # Run UMAP on the patient data
  query <- RunUMAP(query, dims = 1:30, n_neighbors = 20, min_dist = 0.3)

  # Return the processed patient data
  return(query)
}



