#' @title Print PPPR Object Summary
#'
#' @description
#' Make calling the PPPR_OBJECT nicer by printing a summary of its contents.
#'
#' @param x An object of class 'PPPR_OBJECT' to be summarized.
#'
#' @return A summary of the contents of the 'PPPR_OBJECT'.
#' @export
print.PPPR_OBJECT <- function(x) {
  cat("This is an object of class 'PPPR_OBJECT'\n")
  cat("It contains", length(x), "elements:\n\n")

  for (e in seq_along(x)) {
    if (is.list(x[[e]])) {
      cat("Element", e, ":", names(x[e]), "\n")
      cat("A List of", length(x[[e]]),  "matrices\n\n")
    } else if (is.matrix(x[[e]])) {
      cat("Element", e, ":", names(x[e]), "\n")
      cat("A Matrix with dimensions of", paste(dim(x[[e]]), collapse = " x "), "\n\n")
    }
  }
}
