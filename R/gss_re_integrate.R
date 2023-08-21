
#' Re Integrate
#' Splits a Seurat object into a list of Seurat objects to integrate,
#' Performs SCTransform normalization separately for each Seurat object,
#' Runs the PrepSCTIntegration function on the object list,
#' Integrates datasets.
#'
#' @import Seurat
#' @param object Subsetted Seurat object that you want to re-integrate
#' @param ncell_cutoff Minimum number of cells per split object
#'
#' @return Re-integrated Seurat Object.
#' @export
re_integrate <- function(object, ncell_cutoff) {

  ## Split the atlast object by source.
  object.list <- SplitObject(object, split.by = "source")

  ## Loop (backwards) through the `object.list` and remove objects which have less than the chosen cutoff.
  # Start the loop from the last index and iterate backward to the first index
  for (i in length(object.list):1) {
    # Get the number of cells in the Seurat object
    num_cells <- dim(object.list[[i]]@assays$RNA@counts)[2]

    # Check if the number of cells is less than the chosen cutoff.
    if (num_cells < ncell_cutoff) {
      # Remove the Seurat object from the list
      object.list <- object.list[-i]
    }
  }

  ## Run SCTransform on each object in the list
  for (i in 1:length(object.list)) {
    tryCatch({
      object.list[[i]] <- SCTransform(object.list[[i]], method = "glmGamPoi", verbose = FALSE)
    }, error = function(e) {
      message <- paste("Error occurred while applying SCTransform to Seurat object", i, ":")
      message <- paste0(message,
                       e$message,
                       "\n",
                       "To resolve this error Use \`fixInNamespace(\"make_cell_attr\", \"sctransform\")\` to change !identical to !setequal.",
                       "\n",
                       "This is expected to be default in the next update of Seurat",
                       "\n",
                       "\n")
      cat(message)
    })
  }


  ## Select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
  object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
  # Avoid global size limit
  options(future.globals.maxSize = 8000 * 1024^2)
  object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = object.features, verbose = FALSE)

  ## Identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set normalization.method = 'SCT':

  object.anchors <- FindIntegrationAnchors(object.list = object.list, normalization.method = "SCT",
                                                anchor.features = object.features, verbose = FALSE)
  object.integrated <- IntegrateData(anchorset = object.anchors, normalization.method = "SCT",
                                        verbose = FALSE)
  # Now proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset.
  # Commands are identical to the standard workflow, but do not run the ScaleData function after integration.

  object.integrated <- RunPCA(object.integrated, verbose = FALSE)
  object.integrated <- RunUMAP(object.integrated, dims = 1:30)

  ## Plot
  # DimPlot(object.integrated, group.by = c("source", "lv1_annot"), combine = TRUE)


}



