#' pppr_precision
#'
#' @param sample.gs GeneSwitches object that you have already run binarise and GLM on.
#' @param range range of r2cutoffs to use
#'
#' @return the precision dataframe.
#' @export
pppr_precision <- function(sample.gs, range = seq(0.02,0.04,0.005)){

#Build a DF to store accuracy data.
precision <- data.frame(
                r2cutoff = range,                                           # R2cutoff used.
                    n_sg = NA,                                              # Number of genes after min_time_spacing
         inaccuracy_mean = NA                                               # Placeholder for accuracy mean values
)

precision_rownumber <- 1

for (i in range){
  ##Follow the GS and GSS steps.

  # Filter the switching genes
  reference.sg <- filter_switchgenes(sample.gs, allgenes = TRUE, r2cutoff = i)

  # Reduce the binary counts matricies of the query data to only include the selection of evenly distributed genes from the reference.
  sample_reduced      <- filter_gene_expression_for_switching_genes(sample.gs@assays@data@listData$binary   , reference.sg)

  #
  sample.gss <- create_racing_lines(sample_reduced, reference.sg)

  #
  accuracy <- score_gss_accuracy(reference.gss = sample.gss, reference.gs = sample.gs)

  #
  precision$n_sg[precision_rownumber  ] <- dim(reference.sg)[1]
  precision$inaccuracy_mean[precision_rownumber] <- summary(accuracy$inaccuracy)[4]

  #.
  cat(precision_rownumber," ",i," done \n")

  #
  precision_rownumber <- precision_rownumber + 1
}
return(precision)
}


