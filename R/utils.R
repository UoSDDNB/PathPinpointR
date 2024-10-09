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

#' @title welcome_to_PPR
#'
#' @description prints a welcome message when the package is loaded
#'
#' @return "Welcome to PathPinpointR!"
#'
#' @export
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to PathPinpointR!")
}