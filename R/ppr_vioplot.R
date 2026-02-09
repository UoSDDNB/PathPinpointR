#' @title ppr_vioplot
#'
#' @description Creates a violin plot showing the predicted pseudotime
#' positions for each cell in a sample. As well as violin plots of reference cells,
#' subset by ident.
#'
#' @param samples_ppr A list of PPR objects or a single PPR object.
#' @param reference_sce A SingleCellExperiemnt object of the ref trajectory.
#' @param ident selected column from reference_sce colData.
#' 
#' @return Creates a violin plot showing the predicted pseudotime
#' positions for each cell in a sample. As well as violin plots of reference cells,
#' subset by ident.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom vioplot vioplot
#' @importFrom RColorBrewer brewer.pal
#' @export
#'

ppr_vioplot <- function(samples_ppr, reference_sce, ident) {

    ### CHECKS

    # Check if samples_ppr is a single PPR object or a list of them
    if (class(samples_ppr) == "PPR_OBJECT") {
        samples_ppr <- list(samples_ppr)  # Wrap single object in a list
        # name the single object 
        names(samples_ppr) <- "Sample"
    } else if (!all(sapply(samples_ppr, function(x) class(x) == "PPR_OBJECT"))) {
        stop("All elements in samples_ppr must be of class 'PPR_OBJECT'.")
    }

    # check that reference_sce is a SingleCellExperiment object
    if (class(reference_sce) != "SingleCellExperiment") {
        stop("reference_sce must be a SingleCellExperiment object")
    }

    # check that refence_sce has a column named Pseudotime
    if (!"Pseudotime" %in% colnames(colData(reference_sce))) {
        stop("reference_sce must have a column named 'Pseudotime'")
    }

    # check that ident is a column in colData(reference_sce)
    if (!ident %in% colnames(colData(reference_sce))) {
        stop("ident must be a column in colData(reference_sce)")
    }



    ### Data extraction

    # Extract values for violin plot from each PPR object
    vioplot_data <- lapply(samples_ppr, function(ppr) {
        apply(ppr$cells_flat, 1, which_mid_max)
    })
  
    ## Create a data frame to store the pseudotime & idents
    ref_lab_df <- data.frame(
        # Cell names from reference_sce
        cell_names = colnames(reference_sce),
        # True pseudotime values (from reference_sce)
        true_pseudotime = colData(reference_sce)[, "Pseudotime"],
        # True time indices
        true_timeIDX = NA,
        # Reference labels
        idents = droplevels(factor(colData(reference_sce)[, ident]))
    )

    ## Calculate the true time indices based on the pseudotime values
    # Get the pseudotime values from the reference
    ref_ptime <- colData(reference_sce)[, "Pseudotime"]
    # divide the pseudotime range into 100 even intervals
    pt_intervals <- (max(ref_ptime) - min(ref_ptime)) / 100
    # Calculate the true time indices based on the pseudotime values
    ref_lab_df$true_timeIDX <- round((ref_ptime - min(ref_ptime)) / pt_intervals)

    # Loop over each unique value in the 'idents' column
    for (identity in levels(ref_lab_df$idents)) {
        # Filter the rows corresponding to the current 'ident'
        subset_df <- ref_lab_df[ref_lab_df$idents == identity, ]
        
        # Add a new entry to the vioplot_data list
        # The names will be the cell_names, and the values will be the true_timeIDX
        vioplot_data[[identity]] <- setNames(subset_df$true_timeIDX, subset_df$cell_names)
   
    }

    # Calculate mean pseudotime for each sample (including the references)
    mean_pseudotimes <- sapply(vioplot_data, function(data) mean(data, na.rm = TRUE))
    # Order them by their mean pseudotime
    ordered_idents <- names(sort(mean_pseudotimes))
    # Reorder the vioplot_data list based on mean pseudotime
    vioplot_data <- rev(vioplot_data[ordered_idents])


    ### Plotting

    # Set up margins to provide more space for y-axis labels
    # Save current graphical parameters
    old_par <- par(no.readonly = TRUE)
    # Ensure that graphical parameters are reset when the function exits
    on.exit(par(old_par))  # This will restore original settings when the function finishes
    # Adjusting the left margin (2nd value is increased)
    par(mar = c(5, 7.5, 4, 2))  
    
    # Define the colours for the samples
    colours <- brewer.pal(length(vioplot_data), "Paired")
 
    # Define the Y positions for the violin plots
    at_positions <- seq(1, length(vioplot_data), 1)

    # Set up the plot limits and layout
    plot(NULL, xlim = c(0,100), ylim = c(0.5, length(vioplot_data)*1.05),
        xlab = "Pseudotime Index", ylab = "Samples", 
        main = "Pseudotime Distributions", xaxt = "n", yaxt = "n")

    # Add violin plots at specific positions
    for (i in seq_along(vioplot_data)) {
        if (names(vioplot_data)[i] %in% levels(ref_lab_df$idents))
            vioplot(vioplot_data[[i]], horizontal = TRUE, at = at_positions[i], 
            col = "grey", add = TRUE)
        else
        vioplot(vioplot_data[[i]], horizontal = TRUE, at = at_positions[i], 
                col = colours[i], add = TRUE)
    }

    # Add axis labels
    axis(2, at = at_positions, labels = names(vioplot_data), las = 1)

    # add tick to the x-axis
    axis(1, at = seq(0, 100, 10), labels = seq(0, 100, 10))

}

