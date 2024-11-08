#' precisionP
#'
#' @description
#' This function computes the precision (mean inaccuracy) for different values of topnum (number of switching genes)
#' by using a reference sample. It runs in parallel to speed up the computation when large datasets are used.
#' The function returns either a plot of precision (mean inaccuracy) versus the number of switching genes
#' or a data frame containing the computed precision values for each topnum.
#'
#' @param sce A SingleCellExperiment object that has already undergone the following preprocessing steps:
#'   - `GeneSwitches::binarize_exp()`
#'   - `GeneSwitches::find_switch_logistic_fastglm()`
#' 
#' @param n_sg_range A numeric vector specifying the range of `topnum` values to test. Default is `seq(0, 500, 25)`.
#' 
#' @param plot A logical value indicating whether to plot the precision results. Default is `TRUE`. If `FALSE`, the function
#'   will return a data frame with the precision values.
#'
#' @param cpu An integer specifying the number of CPU cores to use for parallel processing. Default is `1`. If `cpu` is 0 or exceeds
#'   the number of available cores, it will default to using all but one core.
#'
#' @return 
#' If `plot = TRUE`, the function returns a plot showing the mean inaccuracy for each value in `n_sg_range`.
#' If `plot = FALSE`, the function returns a data frame with columns `n_sg` (the number of switching genes) and `inaccuracy` (the mean inaccuracy).
#' 
#' @importFrom GeneSwitches filter_switchgenes
#' @importFrom graphics points grid legend
#' @importFrom methods is
#' @importFrom parallel makeCluster detectCores clusterExport parLapply stopCluster
#'
#' @export
precisionP <- function(sce,
                      n_sg_range = seq(0, 500, 25),
                      plot = TRUE,
                      cpu = 1) {

  # Check if sce is a SingleCellExperiment object
  if (!is(sce, "SingleCellExperiment")) {
    stop(paste0("\"", deparse(substitute(sce)), "\" must be a SingleCellExperiment object."))
  }

  # Check if sce contains binarized expression data
  if (is.null(sce@assays@data@listData$binary)) {
    stop("\n  Binarized expression data not found in \"",
         deparse(substitute(sce)),
         "\" \n  You must run `GeneSwitches::binarize_exp()` first.")
  }

  # Set cpu to all but one core if cpu is 0 or exceeds available cores
  if (cpu == 0 || cpu >= detectCores()) {
    cpu <- detectCores() - 1
  }

  # Set up a parallel cluster
  library(parallel)
  cl <- makeCluster(cpu)

  # Export necessary objects and functions to each worker
  clusterExport(cl,
                c("sce",
                  "filter_switchgenes",
                  "predict_position",
                  "accuracy_test",
                  "which_mid_max",
                  "rowData"),  # Added which_mid_max here
                envir = environment())

  # Parallel execution using parLapply
  precision_results <- parLapply(cl, n_sg_range, function(i) {
    # Filter switching genes and perform position prediction and accuracy testing
    switching_genes <- filter_switchgenes(sce, allgenes = TRUE, r2cutoff = 0, topnum = i)
    sample_ppr <- predict_position(sce, switching_genes)
    accuracy <- accuracy_test(ppr = sample_ppr, reference_sce = sce, plot = FALSE)

    # Collect results for each iteration
    data.frame(n_sg = i,
               inaccuracy = mean(accuracy$inaccuracy))
  })

  # Stop the cluster after all tasks are completed
  stopCluster(cl)

  # Combine results into a final data frame
  precision <- do.call(rbind, precision_results)
  precision$n_sg <- n_sg_range

  # Plot if requested
  if (plot) {
    plot(precision$n_sg, precision$inaccuracy, type = "l",
         xlab = "Number of Switching Genes", ylab = "Mean Inaccuracy",
         main = "Mean Inaccuracy by Number of Switching Genes")
    points(precision$n_sg, precision$inaccuracy, pch = 16)
    text(x = precision$n_sg[which.min(precision$inaccuracy)],
         y = precision$inaccuracy[which.min(precision$inaccuracy)],
         labels = precision$n_sg[which.min(precision$inaccuracy)], pos = 1)
    grid()
    myplot <- recordPlot()
    return(myplot)
  } else {
    return(precision)
  }
}