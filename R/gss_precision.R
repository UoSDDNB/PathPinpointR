#' gss_precision
#'
#' @param sample.gs GeneSwitches object that you have already run binarise and GLM on.
#' @param range range of r2cutoffs to use
#'
#' @return the precision dataframe.
#' @export
gss_precision <- function(sample.gs, range = seq(0.02,0.04,0.005)){

#Build a DF to store accuracy data.
precision2 <- data.frame(
  r2cutoff = range,
  n_sg          = NA,                                                   # Number of genes after min_time_spacing
  inaccuracy_min     = NA,                                              # Placeholder for accuracy minimum values
  inaccuracy_1       = NA,
  inaccuracy_median  = NA,
  inaccuracy_mean    = NA,
  inaccuracy_3       = NA,
  inaccuracy_max     = NA
)

precision_rownumber <- 1

for (i in range){

  reference.sg <- filter_switchgenes(sample.gs, allgenes = TRUE, r2cutoff = i)

  # Reduce the binary counts matricies of the query data to only include the selection of evenly distributed genes from the refernence.
  sample_reduced      <- filter_gene_expression_for_switching_genes(sample.gs@assays@data@listData$binary   , reference.sg)

  #
  sample.gss <- create_racing_lines(sample_reduced, reference.sg)

  #
  accuracy <- score_gss_accuracy(reference.gss = sample.gss, reference.gs = sample.gs)

  #
  precision2$n_sg[precision_rownumber  ] <- dim(reference.sg)[1]
  precision2$inaccuracy_min[precision_rownumber ] <- summary(accuracy$inaccuracy)[1]
  precision2$inaccuracy_1[precision_rownumber ] <- summary(accuracy$inaccuracy)[2]
  precision2$inaccuracy_median[precision_rownumber] <- summary(accuracy$inaccuracy)[3]
  precision2$inaccuracy_mean[precision_rownumber] <- summary(accuracy$inaccuracy)[4]
  precision2$inaccuracy_3[precision_rownumber ] <- summary(accuracy$inaccuracy)[5]
  precision2$inaccuracy_max[precision_rownumber ] <- summary(accuracy$inaccuracy)[6]

  precision_rownumber <- precision_rownumber + 1
  #
  cat(i," done \n")
}

# ggplot(data = precision2, mapping = aes(x = r2cutoff)) +
#   geom_point(aes(y = inaccuracy_mean), color = "blue")
#
# ggplot(data = precision2, mapping = aes(x = n_sg)) +
#   geom_point(aes(y = inaccuracy_mean), color = "blue") +
#   scale_x_log10()
#
# ggplot(data = precision2, mapping = aes(x = min_time_spacing)) +
#   geom_point(aes(y = n_genes), color = "red") +
#   scale_y_log10()


return(precision2)
}


