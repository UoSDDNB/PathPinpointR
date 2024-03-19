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
      cat("Element", e, ":", names(x[e]), "\n")
      cat("A List of", length(x[[e]]),  "matrices\n\n")
    } else if (is.matrix(x[[e]])) {
      cat("Element", e, ":", names(x[e]), "\n")
      cat("A Matrix with dimensions of", paste(dim(x[[e]]), collapse = " x "), "\n\n")
    } else {
      cat("Element", e, ":", names(x[e]), "\n")
      cat("Type:", class(x[[e]]), "\n\n")
    }
  }
}

#' Set Default Print Method for PPR_OBJECT
setMethod("print", "PPR_OBJECT", print.PPR_OBJECT)
