# R/utils.R

#' PPR_OBJECT Class Definition
#'
#' This class is used within PathPinpointR to store and manage predicted positions.
#' @export
setClass("PPR_OBJECT")

#' @title Print method for PPR Object.
#'
#' @description
#' Make calling the PPR_OBJECT nicer by printing a summary of its contents.
#' This function can be used to quickly inspect the structure of your results.
#'
#' @param x An object of class 'PPR_OBJECT' to be summarized.
#'
#' @return A summary of the contents of the object.
#' @export
print.PPR_OBJECT <- function(x) {
  if (!inherits(x, "PPR_OBJECT")) {
    stop("Input is not of class 'PPR_OBJECT'")
  }

  cat("This is an object of class 'PPR_OBJECT'\n")
  cat("It contains", length(x), "elements:\n\n")

  for (e in seq_along(x)) {
    if (is.list(x[[e]])) {
      cat("Element",
          e,
          ":",
          names(x[e]),
          "\n")
      cat(" A List of",
          length(x[[e]]),
          "matrices\n")
      cat(" Each of dimensions",
          paste(dim(x[[e]][[1]]),
          collapse = " x "),
          "\n\n")
    } else if (is.matrix(x[[e]])) {
      cat("Element",
          e,
          ":",
          names(x[e]),
          "\n")
      cat(" A Matrix with dimensions of",
          paste(dim(x[[e]]),
          collapse = " x "),
          "\n\n")
    } else if (names(x)[e] == "sd") {
      cat("Standard Deviation = ",
          x[[e]],
          "\n")
    } else if (names(x)[e] == "z_score") {
      cat("Z-Score            = ",
          x[[e]],
          "\n")
    } else if (names(x)[e] == "p_value") {
      cat("p-Value            = ",
          x[[e]],
          "\n\n")
    } else {
      cat("Element",
          e,
          ":",
          names(x[e]),
          "\n")
      cat("Type:",
          class(x[[e]]),
          "\n\n")
    }
  }
}

#' Set Default Print Method for PPR_OBJECT
setMethod("print", "PPR_OBJECT", print.PPR_OBJECT)


#' @title Select the mid index among multiple occurrences of the max value.
#'
#' @description This helper function is used within PathPinpointR to resolve ties when
#' multiple values are equal to the maximum.
#'
#' @param n Numeric vector.
#' @return The middle index among multiple occurrences of the maximum value.
#' @export
which_mid_max <- function(n) {
  max_indices <- which(n == max(n))
  middle_index <- ceiling(length(max_indices) / 2)
  max_indices[middle_index]
}


#' @title get_example_data
#'
#' @description Downloads example data from Dropbox. 
#'
#' @return Example data saved to the current working directory.
#' @export
get_example_data <- function() {
  dest_files <- c("./reference.rds", "./LW120.rds", "./LW122.rds")
  dropbox_urls <- c("https://www.dropbox.com/scl/fi/bfklckyfsxv6tswejhk9u/blastocyst_downsampled.rds?rlkey=lwnccfuhf8beq3viupg1l41lv&st=w9ncvlt9&dl=1", # nolint: line_length_linter.
                    "https://www.dropbox.com/scl/fi/eo397e2whn4gb6wasn66l/LW120_seu.rds?rlkey=hr2ihl8iul4tye1tccznb62yn&st=00jmtbzz&dl=1", # nolint: line_length_linter.
                    "https://www.dropbox.com/scl/fi/088eocoyn5xrxc1ngvo4d/LW122_seu.rds?rlkey=cxfhlg835324gl00bh4yidl57&st=t5vdapbj&dl=1") # nolint: line_length_linter.

  for (i in seq_along(dest_files)) {
    if (!file.exists(dest_files[i])) {
      download_command <- paste("curl -L -o",
                                shQuote(dest_files[i]),
                                shQuote(dropbox_urls[i]))

      exit_status <- system(download_command)

      if (exit_status == 0) {
        print(paste("File successfully downloaded to", dest_files[i]))
      } else {
        print("Error: Failed to download the file.")
      }
    } else {
      print(paste("File already exists at", dest_files[i]))
    }
  }
}

#' @title get_synthetic_data
#'
#' @description Generates a synthetic trajectory and two subsets of the trajectory
#' to be used as example samples.
#'
#' @return A list of three SingleCellExperiment objects: one reference and two samples.
#' 
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @export
get_synthetic_data <- function() {
  # number of cells 
  n_cells <- 300

  # generate the means of the gene expression
  means <- rbind(
      # non-DE genes
      matrix(rep(rep(c(0.1,0.5,1,2,3), each = n_cells),100),
          ncol = n_cells, byrow = TRUE),
      # deactivation 1
      matrix(rep((exp(atan((n_cells:1) - 295) * 4)), 150), ncol = n_cells, byrow = TRUE),
      # deactivation 2
      matrix(rep(exp(atan(((n_cells:1) - 290) * 2)), 100), ncol = n_cells, byrow = TRUE),
      # deactivation 3
      matrix(rep(exp(atan(((n_cells:1) - 280) / 2)), 50), ncol = n_cells, byrow = TRUE),
      # deactivation 4
      matrix(rep(exp(atan(((n_cells:1) - 250) / 10)), 50), ncol = n_cells, byrow = TRUE),
      # deactivation 5
      matrix(rep(exp(atan(((n_cells:1) - 200) / 50)), 50), ncol = n_cells, byrow = TRUE),
      # deactivation 6
      matrix(rep(exp(atan(((n_cells:1) - 150) / 50)), 50), ncol = n_cells, byrow = TRUE),
      # deactivation 7
      matrix(rep(exp(atan(((n_cells:1) - 100) / 50)), 50), ncol = n_cells, byrow = TRUE),
      # deactivation 8
      matrix(rep(exp(atan(((n_cells:1) - 50) / 50)), 50), ncol = n_cells, byrow = TRUE),
      # deactivation 9
      matrix(rep(exp(atan(((n_cells:1) - 25) / 50)), 50), ncol = n_cells, byrow = TRUE),
      # deactivation 10
      matrix(rep(exp(atan(((n_cells:1)) - 0 / 50)), 50), ncol = n_cells, byrow = TRUE),

      # activation 1
      matrix(rep(exp(atan(((1:n_cells) - 10) / 50 )), 50), ncol = n_cells, byrow = TRUE),
      # activation 2
      matrix(rep(exp(atan(((1:n_cells) - 50) / 50 )), 50), ncol = n_cells, byrow = TRUE),
      # activation 3
      matrix(rep(exp(atan(((1:n_cells) - 100) / 50 )), 50), ncol = n_cells, byrow = TRUE),
      # activation 4
      matrix(rep(exp(atan(((1:n_cells) - 150) / 50 )), 50), ncol = n_cells, byrow = TRUE),
      # activation 5
      matrix(rep(exp(atan(((1:n_cells) - 200) / 50 )), 50), ncol = n_cells, byrow = TRUE),
      # activation 6
      matrix(rep(exp(atan(((1:n_cells) - 250) / 50 )), 50), ncol = n_cells, byrow = TRUE),
      # activation 7
      matrix(rep(exp(atan(((1:n_cells) - 270) / 10 )), 50), ncol = n_cells, byrow = TRUE),
      # activation 8
      matrix(rep(exp(atan(((1:n_cells) - 280) / 2 )), 50), ncol = n_cells, byrow = TRUE),
      # activation 9
      matrix(rep(exp(atan(((1:n_cells) - 290) * 2  )), 100), ncol = n_cells, byrow = TRUE),
      # activation 10
      matrix(rep(exp(atan(((1:n_cells) - 295) * 4  )), 150), ncol = n_cells, byrow = TRUE),

      # transient
      matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),55), 
        ncol = n_cells, byrow = TRUE)
    )

  # set seed for reproducaibility
  set.seed(2)

  # simulate a counts matrix
  counts <- apply(means,2,function(cell_means){
      # total number of counts
      total <- rnbinom(1, mu = 7500, size = 4)
      # simulate counts
      rmultinom(1, total, cell_means)
  })

  # name the genes and cells G and c respectively
  rownames(counts) <- paste0('G',1:dim(means)[1] )
  colnames(counts) <- paste0('c',1:n_cells )

  # assign the counts to a SingleCellExperiment object
  reference_sce <- SingleCellExperiment(assays = List(counts = counts))

  # filter out lowly expressed genes 
  geneFilter <- apply(assays(reference_sce)$counts,1,function(x){
      sum(x >= 3) >= 10
  })
  reference_sce <- reference_sce[geneFilter, ]

  # normalize the data
  FQnorm <- function(counts){
      rk <- apply(counts,2,rank,ties.method='min')
      counts.sort <- apply(counts,2,sort)
      refdist <- apply(counts.sort,1,median)
      norm <- apply(rk,2,function(r){ refdist[r] })
      rownames(norm) <- rownames(counts)
      return(norm)
  }
  assays(reference_sce)$norm <- FQnorm(assays(reference_sce)$counts)
  # rename for compatability with GS
  assays(reference_sce)$expdata <- assays(reference_sce)$norm

  # PCA
  pca <- prcomp(t(log1p(assays(reference_sce)$norm)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  reducedDims(reference_sce) <- SimpleList(PCA = rd1)

  # clustering 
  cl1 <- Mclust(rd1)$classification
  colData(reference_sce)$GMM <- cl1

  # name the clusters
  cluster_names <- c("1",
                    "2",
                    "3",
                    "4",
                    "5",
                    "6")
                    
  colData(reference_sce)$clust_names <- factor(cluster_names[colData(reference_sce)$GMM],
                                              levels = cluster_names)

  # return the reference sce
  return(reference_sce)
}

#' @title welcome_to_PPR
#'
#' @description prints a welcome message when the package is loaded
#' 
#' @return "Welcome to PathPinpointR!"
#' @export
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to PathPinpointR!")
}

#' @title ppr2seu
#'
#' @description This wrapper function extracts the predicted position of each 
#' cell from a ppr object and adds it to the metdata of a Seurat object.
#'
#' @param ppr A ppr object.
#' @param seu A Seurat object.
#' @param colname The name of the column to add to the Seurat object.
#' 
#' @return The Seurat object with the predicted position added to the metadata.
#' 
#' @importFrom Seurat AddMetaData
#' 
#' @export
ppr2seu <- function(ppr, seu, colname = "PPR_pseudotime_idx") {
    # use apply and which_mid_max to find the index of the max value in each row
    predicted_ptime <- apply(ppr$cells_flat, 1, which_mid_max)
    # predicted_ptime is an integer with names corresponding to the cell names

    # Ensure that all cell names in predicted_ptime exist in the Seurat object
        if (!identical(names(predicted_ptime), colnames(seu))) {
            stop("Mismatch between cell names in PPR and Seurat object")
        }

    # add the predicted positions to the metadata of the seurat object
    seu <- AddMetaData(seu, predicted_ptime, colname)

    return(seu)
}