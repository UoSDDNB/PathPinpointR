#' precision
#'
#' @description
#' Used to find the optimum r2cutoff to use,
#' #' based on using reference as your sample.
#'
#' @param sce SingleCellExperiement object,
#'  that you have already run GeneSwitches::binarize_exp(),
#'  and GeneSwitches::find_switch_logistic_fastglm() on.
#' @param r2_cutoff_range range of r2cutoffs to use
#' @param plot logical, do you want to plot the precision df.
#'
#' @return a plot of the precision df.
#' @importFrom GeneSwitches filter_switchgenes
#' @importFrom graphics points grid legend
#' @importFrom methods is
#'
#' @export
precision <- function(sce,
                      r2_cutoff_range = seq(0.0, 0.5, 0.1),
                      plot = "r2") {

  # Check that a suitable plot has been selected
  if (!plot %in% c("r2", "n_sg", FALSE)) {
    stop("plot must be either 'r2' or 'n_sg'")
  }

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
    r2cutoff = r2_cutoff_range,      # R2cutoff used.
    n_sg = NA,                   # Number of switching genes
    inaccuracy = NA
  )

  nrow_precision <- 1

  for (i in r2_cutoff_range){
    ##Follow the GS and PPR steps.

    # Filter the switching genes
    switching_genes <- filter_switchgenes(sce, allgenes = TRUE, r2cutoff = i)

    # Reduce the binary counts matricies of the query data,
    # to only include the selection of switching genes from the reference.
    sample_reduced <- subset_switching_genes(sce, switching_genes)

    #
    sample_ppr <- predict_position(sample_reduced, switching_genes)

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



  if (plot == "r2" || plot == "n_sg") {
    if (plot == "r2") {
      # Plot r2 by inaccuracy
      plot(precision$r2cutoff,
           precision$inaccuracy,
           type = "l",
           xlab = "R-squared Cutoff",
           ylab = "Mean Inaccuracy",
           main = "Mean Inaccuracy vs R-squared Cutoff")
      # Adding points
      points(precision$r2cutoff, precision$inaccuracy, pch = 16)
    } else if (plot == "n_sg") {
      # Plot n_sg by inaccuracy
      plot(precision$n_sg, precision$inaccuracy, type = "l",
           xlab = "Number of Switching Genes", ylab = "Mean Inaccuracy",
           main = "Mean Inaccuracy by Number of Switching Genes")
      # Adding points
      points(precision$n_sg, precision$inaccuracy, pch = 16)

      # label the point at the lowest inaccuracy
      text(x = precision$n_sg[which.min(precision$inaccuracy)],
           y = precision$inaccuracy[which.min(precision$inaccuracy)],
           labels = paste("R^2 =",
                          precision$r2cutoff[which.min(precision$inaccuracy)]),
           pos = 1)
    }

    # Adding grid
    grid()

    # Save the plot as an object
    myplot <- recordPlot()
    return(myplot)
  } else {
    return(precision)
  }

}