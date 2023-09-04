#' #' gss_plot_cell_range
#' #'
#' #' @param sample.gss GSS OBJECT
#' #' @param cell_idxs Range or list of cell indexes
#' #' @param genes_of_interest The genes that you would like to plot
#' #' @param switching_genes The data which includes all of the switching genes
#' #'
#' #' @return plot to show a large number of cell positions
#' #' @export
#' gss_plot_cell_range <- function(sample.gss, cell_idxs = c(1:2,4), genes_of_interest = NULL, switching_genes = NULL) {
#'   for (i in cell_idxs){
#'     if (i > 1) {
#'       lines(x = 1:100,
#'             y = sample.gss$cells_flat[i,],
#'             ylim = c(0,max(sample.gss$cells_flat[i,]) + 10),
#'             xlim = c(0,100),
#'             pch = 20,
#'             cex = 0.8,
#'             type = "l")
#'     } else if (i == 1) {
#'       plot(x = 1:100,
#'            y = sample.gss$cells_flat[i,],
#'            ylim = c(0,max(sample.gss$cells_flat[i,]) + 10),
#'            xlim = c(0,100),
#'            pch = 20,
#'            cex = 0.8,
#'            type = "l",
#'            xlab = "Pseudo-Time Index",
#'            ylab = "GSS Score",
#'            main = paste("Cell Positions"))
#'
#'         }
#'   }
#'
#'   if (length(genes_of_interest) > 0) {
#'     for (gene_name in genes_of_interest) {
#'       segments(switching_genes[gene_name,"switch_at_timeidx"],
#'                -3.9,
#'                switching_genes[gene_name,"switch_at_timeidx"],
#'                -0.5,
#'                lwd = 1,
#'                lty = 2)
#'
#'       text(x = switching_genes[gene_name,"switch_at_timeidx"] + 3,
#'            y = 1,
#'            labels = gene_name,
#'            srt = -20,
#'            cex = 0.86,
#'            pos = 2)
#'     }
#'   }
#' }
