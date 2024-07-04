# R/utils.R

#' PPR_OBJECT Class Definition
#'
#' @export
setClass("PPR_OBJECT")

#' @title Print method for PPR Object.
#'
#' @description
#' Make calling the PPR_OBJECT nicer by printing a summary of its contents.
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


#' @title Select the middle index among multiple occurrences of the maximum value
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
#' @description downloads example data from dropbox
#'
#' @return example data is saved to the current working directory
#'
#' @export
#'
get_example_data <- function() {
  dest_file <- "./reference.rds"
  if (!file.exists(dest_file)) {
    dropbox_url <- "https://www.dropbox.com/scl/fi/t9gaxkxb97adoxgemou9o/binarized_Petro16_sce.rds?rlkey=f3rm20bk0s41138wq7422y55p&dl=1"
    download_command <- paste("curl -L -o",
                              shQuote(dest_file),
                              shQuote(dropbox_url))

    exit_status <- system(download_command)

    if (exit_status == 0) {
      print(paste("File successfully downloaded to", dest_file))
    } else {
      print("Error: Failed to download the file.")
    }
  } else {
    print(paste("File already exists at", dest_file))
  }
}
