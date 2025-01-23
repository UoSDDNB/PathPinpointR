#' precision
#'
#' @description
#' Used to find the optimum number of switching genes. This function is often
#' iterated over in the PathPinpointR pipeline. Refer to "README.Rmd" for an
#' example workflow using different topnum values.
#'
#' @param sce SingleCellExperiement object,
#'  that you have already run GeneSwitches::binarize_exp(),
#'  and GeneSwitches::find_switch_logistic_fastglm() on.
#' @param n_sg_range range of topnum values to use
#' @param plot logical, do you want to plot the precision df.
#'
#' @return a plot of the precision df.
#' @importFrom GeneSwitches filter_switchgenes
#' @importFrom graphics points grid legend
#' @importFrom methods is
#'
#' @export
precision <- function(sce,
                      n_sg_range = seq(0, 500, 25),
                      plot = TRUE) {

  # Check if sce is a SingleCellExperiment object
  if (!is(sce, "SingleCellExperiment")) {
    # If sce is not a SingleCellExperiment object, display a message and stop.
    stop(paste0("\"",
                deparse(substitute(sce)),
                "\" must be a SingleCellExperiment object."))
  }

  # Check if sce contains binarized expression data
  if (is.null(sce@assays@data@listData$binary)) {
    # If binarized expression data is missing, display a message and stop.
    stop("\n  Binarized expression data not found in \"",
         deparse(substitute(sce)),
         "\" \n  You must run `GeneSwitches::binarize_exp()` first.")
  }

  # Build a DF to store accuracy data.
  precision <- data.frame(
    n_sg = n_sg_range,   # number of switching genes
    inaccuracy = NA
  )

  nrow_precision <- 1

  for (i in n_sg_range){
    ##Follow the GS and PPR steps.

    # Filter the switching genes
    switching_genes <- filter_switchgenes(sce,
                                          allgenes = TRUE,
                                          r2cutoff = 0,
                                          topnum = i)

    #
    sample_ppr <- predict_position(sce, switching_genes)

    #
    accuracy <- accuracy_test(ppr = sample_ppr,
                              reference_sce = sce,
                              plot = FALSE)

    #
    precision$n_sg[nrow_precision] <- dim(switching_genes)[1]
    precision$inaccuracy[nrow_precision] <- summary(accuracy$inaccuracy)[4]

    #
    cat(nrow_precision, " ", i, " done \n")

    #
    nrow_precision <- nrow_precision + 1
  }



  if (plot) {
    # Plot n_sg by inaccuracy
    plot(precision$n_sg, precision$inaccuracy, type = "l",
         xlab = "Number of Switching Genes", ylab = "Mean Inaccuracy",
         main = "Mean Inaccuracy by Number of Switching Genes")
    # Adding points
    points(precision$n_sg, precision$inaccuracy, pch = 16)

    # label the point at the lowest inaccuracy
    text(x = precision$n_sg[which.min(precision$inaccuracy)],
         y = precision$inaccuracy[which.min(precision$inaccuracy)],
         labels = precision$n_sg[which.min(precision$inaccuracy)],
         pos = 1)

    # Adding grid
    grid()

    # Save the plot as an object
    myplot <- recordPlot()
    return(myplot)
  } else {
    return(precision)
  }

}
