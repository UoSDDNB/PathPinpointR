#' #' @title Select Evenly Distributed Switching Genes
#' #'
#' #' @description
#' #' Filters switching genes from GenesSwitches.
#' #'
#' #' @param sg_allgenes Filtered Switching Genes from Gene Switches. (dont use `topnum` instead choose a Psvalue cutoff )
#' #' @param min_time_spacing The minimum distance between genes when evenly distributing.
#' #'
#' #' @return Switching genes which are evenly distributed through pseudotime at a selected density.
#' #' @export
#' select_evenly_distributed_switching_genes <- function(sg_allgenes, min_time_spacing = 1){
#'
#'   ## Checks
#'   # Make sure there is a minimum of 1 gene switched on (up)
#'   stopifnot("Too few genes switching \"up\" " = !sum(sg_allgenes$direction == "up") < 2)
#'   # Make sure there is a minimum of 1 gene switched off (down)
#'   stopifnot("Too few genes switching \"down\" " = !sum(sg_allgenes$direction == "down") < 2)
#'
#'
#'   ## Sort sg_allgenes by pseudoR2s
#'   sg_allgenes <- sg_allgenes[order(-sg_allgenes$pseudoR2s),]
#'
#'   # Initialize the time value of the previously selected gene for both "up"'s and "down"'s
#'   ups <- sg_allgenes[sg_allgenes$direction == "up", ]
#'   prev_ups_times <- ups$switch_at_timeidx[1]
#'   downs <- sg_allgenes[sg_allgenes$direction == "down", ]
#'   prev_downs_times <- downs$switch_at_timeidx[1]
#'
#'   # Initialize the subsetted matrix with the first "up" and "down" gene
#'   subsetted_matrix <- downs[1, ]
#'   subsetted_matrix <- rbind(subsetted_matrix, ups[1, ])
#'
#'   # Loop over the remaining "up" genes and add them to the subsetted matrix if they meet the criteria
#'   for (i in 2:nrow(ups)) {
#'     # Check if the time value of the current gene is spaced by at least min_time_spacing from all previously selected genes
#'     if (all(abs(ups$switch_at_timeidx[i] - prev_ups_times) >= min_time_spacing)) {
#'       # Add the current gene to the subsetted matrix
#'       subsetted_matrix <- rbind(subsetted_matrix, ups[i, ])
#'       # Update the previous time values
#'       prev_ups_times <- c(prev_ups_times, ups$switch_at_timeidx[i])
#'     }
#'   }
#'
#'   # Loop over the remaining "down" genes and add them to the subsetted matrix if they meet the criteria
#'   for (i in 2:nrow(downs)) {
#'     # Check if the time value of the current gene is spaced by at least min_time_spacing from all previously selected genes
#'     if (all(abs(downs$switch_at_timeidx[i] - prev_downs_times) >= min_time_spacing)) {
#'       # Add the current gene to the subsetted matrix
#'       subsetted_matrix <- rbind(subsetted_matrix, downs[i, ])
#'       # Update the previous time values
#'       prev_downs_times <- c(prev_downs_times, downs$switch_at_timeidx[i])
#'     }
#'   }
#'
#'   # return the subsetted matrix
#'   gs_scorer_genes <- subsetted_matrix
#'   return(gs_scorer_genes)
#' }
