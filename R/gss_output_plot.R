#' GSS OUTPUT PLOT
#'
#' @description
#' Plots the predicted positon of your sample.
#'
#' @param sample.gss A GSS_OBJECT of the sample you wish to plot
#' @param col The colour that you'd like
#' @param overlay set to TRUE if you would like this plot to overlay a previous plot.
#'
#' @return nice plot highlighting the probable position of your sample on your trajectory.
#' @export
#'
gss_output_plot <- function(sample.gss, col = "red", overlay = FALSE){

 if (!overlay) {
  plot(x = 1:100,
       y = (sample.gss$sample_flat/max(sample.gss$sample_flat)*100),
       ylim = c(0,100),
       pch = 20,
       cex = 0.8,
       col = col,
       type = "p",
       xlab = "Pseudo-Time Index",
       ylab = "GSS Score",
       main = "Predicted Position of Sample")

  segments(which.max(colSums(sample.gss$sample_flat)),
           -3.9,
           which.max(colSums(sample.gss$sample_flat)),
           (sample.gss$sample_flat[which.max(colSums(sample.gss$sample_flat))]/max(sample.gss$sample_flat)*100),
           lwd = 1,
           lty = 2,
           col = col)
 } else {
   lines( x = 1:100,
          y = (sample.gss$sample_flat/max(sample.gss$sample_flat)*100),
        pch = 20,
        cex = 0.8,
        col = col,
        type = "p")

   segments(which.max(colSums(sample.gss$sample_flat)),
            -3.9,
            which.max(colSums(sample.gss$sample_flat)),
            (sample.gss$sample_flat[which.max(colSums(sample.gss$sample_flat))]/max(sample.gss$sample_flat)*100),
            lwd=1,
            lty = 2,
            col = col)
   }

}



