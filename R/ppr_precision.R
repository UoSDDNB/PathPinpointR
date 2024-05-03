#' ppr_precision
#'
#' @description Used to find the optimum r2cutoff to use, based on using reference as your sample.
#' 
#' @param sample.gs GeneSwitches object that you have already run binarise and GLM on.
#' @param r2_cutoff_range range of r2cutoffs to use
#'
#' @return a plot of the precision df. (#TODO make the df as an optional output!)
#' @importFrom GeneSwitches filter_switchgenes
#' @export
ppr_precision <- function(sample.gs, r2_cutoff_range = seq(0.0,0.5,0.1)){

#Build a DF to store accuracy data.
precision <- data.frame(
                r2cutoff = r2_cutoff_range,                                           # R2cutoff used.
                    n_sg = NA,                                              # Number of genes after min_time_spacing
         inaccuracy_mean = NA                                               # Placeholder for accuracy mean values
)

precision_rownumber <- 1

for (i in r2_cutoff_range){
  ##Follow the GS and GSS steps.

  # Filter the switching genes
  reference.sg <- filter_switchgenes(sample.gs, allgenes = TRUE, r2cutoff = i)

  # Reduce the binary counts matricies of the query data to only include the selection of evenly distributed genes from the reference.
  sample_reduced      <- ppr_filter_gene_expression_for_switching_genes(sample.gs@assays@data@listData$binary   , reference.sg)

  #
  sample.ppr <- ppr_predict_position(sample_reduced, reference.sg)

  #
  accuracy <- ppr_accuracy_test(reference.ppr = sample.ppr, reference.gs = sample.gs, plot = FALSE)

  #
  precision$n_sg[precision_rownumber  ] <- dim(reference.sg)[1]
  precision$inaccuracy_mean[precision_rownumber] <- summary(accuracy$inaccuracy)[4]

  #.
  cat(precision_rownumber," ",i," done \n")

  #
  precision_rownumber <- precision_rownumber + 1
}

# Plotting inaccuracy_mean
plot(precision$r2cutoff, precision$inaccuracy_mean, type = "l",  # Line plot
     xlab = "R-squared Cutoff", ylab = "Mean Inaccuracy",        # Axis labels
     main = "Mean Inaccuracy vs R-squared Cutoff")               # Plot title

# Adding points
points(precision$r2cutoff, precision$inaccuracy_mean, pch = 16)

# Adding grid
grid()

# Adding a legend
legend("topright", legend = "Inaccuracy Mean", pch = 16, col = "black", bty = "n")

# Save the plot as an object
myplot <- recordPlot()


return(myplot)
}


