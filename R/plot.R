#' Visualize gene expression profile of reference
#'
#' @description Given current cell type x gene reference matrix, visualize it
#' using heatmap. Allow visualization to be done in log scale. In default,
#' it automatically reorders genes for better interpretation.
#'
#' @param ref reference matrix
#' @param want_log whether visualization should be done in log scale (Default: F)
#' @param cluster_col whether genes should be reordered (Default: T)
#' @param cluster_row whether cell types should be reordered (Default: F)
#' @param title Default: Gene Expression Heatmap
#'
#' @return Heatmap visualisation of gene expression profiles of reference
#'
#' @examples
#'
#' data(mouse_hypothalamus)
#' new_reference <- appendPlaceholder(mouse_hypothalamus$reference, numPlaceholder=1)
#' p <- visualizeReference(new_reference)
#' p
#'
#' @export


visualizeReference <- function(ref, want_log=F, cluster_col=T, cluster_row=F,
                                title = "Gene Expression Heatmap") {


  if (want_log) {
    new_ref = log(ref)
  } else {
    new_ref = ref
  }

  min_val <- min(new_ref)
  max_val <- max(new_ref)

  heatmap_plot <- pheatmap::pheatmap(new_ref,
                                     cluster_rows = F,
                                     cluster_cols = cluster_col,
                                     treeheight_row = 0,
                                     treeheight_col = 0,
                                     show_rownames = T,
                                     show_colnames = F,
                                     color = colorRampPalette(c("white", "orange", "red"))(50),
                                     main = title,
                                     breaks = seq(min_val, max_val, length.out = 51))


  return(heatmap_plot)
}





#' Visualizing deconvolution output for Spatial Transcriptomics data that have
#' multiple slices
#'
#' @description Given the deconvolution output for Spatial Transcriptomics, plot
#' the result to a grid, every entry represents a single slice
#'
#' @param deconvResult deconvolution result, dimension: pixel x cell type
#' @param meta metadata for the Spatial Transcriptomics data, row names must be
#'    barcodes, it also should have x, y as two columns, as well as a column for slice
#' @param color named vector representing the color you want for each cell type
#' @param numRow number of row of grid
#' @param numCol number of column of grid
#' @param r radius for each pixel as a circle (Default: 40)
#' @param title title of the graph (Default: "Cell Type Proportion per Pixel")
#' @param legend_grob object that will be displayed as legend
#' @param verbose Default: T
#' @param topic set a number if you want to visualize a specific cell type (Default: NA)
#'
#' @return Piechart visualization for deconvolution result
#'
#' @examples
#'
#' data(mouse_hypothalamus)
#' color = c('Astrocyte'='red', 'Endothelial'='orange',
#' 'Ependymal'='black', 'Excitatory'='blue',
#' 'Inhibitory'='green', 'Microglia'='purple',
#' 'OD Immature'='yellow', 'OD Mature'='brown', 'Pericytes'='cyan')
#' legend_labels <- names(color)
#' legend_colors <- as.vector(color)
#' legend_grob <- getLegendGrob(legend_labels, legend_colors)
#' plotDeconvRes(mouse_hypothalamus$ground_truth_deconv,
#' mouse_hypothalamus$spatial_meta,
#' color = color,
#' numRow = 1,
#' numCol = 2,
#' legend_grob = legend_grob)
#'
#' @export

plotDeconvRes <- function(deconvResult, meta, color, numRow, numCol, r=40,
                          title = "Cell Type Proportion per Pixel",
                          legend_grob=NA, verbose=T, topic=NA) {
  ordered_col = color[order(names(color))]
  plt_lst <- plot_est_pixels(deconvResult, meta, ordered_col, verbose, topic,r=r)
  plot_a_grid(plt_lst, legend_grob = legend_grob, num_row = numRow,
              num_col = numCol, title=title)
}





#' Helper function to get a list of deconvolution results, one for each slice
#'
#' @description For each slice, we want to plot a pie chart and collect it
#' within a list
#'
#' @param deconvResult deconvolution result for all slices, in matrix form with
#' dimension: pixel x cell type
#' @param meta metadata for Spatial Trasnscriptomics data, require slice column,
#' x and y coordinates and pixel ID column (barcodes)
#' @param color color vector named by different cell types
#' @param verbose Default: T
#' @param topic if you want to visualize a specific cell type, enter a number
#' for this variable
#' @param r radius for each pixel (Default: 40)
#'
#' @return Pie charts for each slice in a list
#'
#' @noRd
#'
#' @keywords internal

plot_est_pixels <- function(deconvResult, meta, color, verbose=T, topic=NA, r=40) {

  plot_list = list()
  for (slice in unique(meta$slice)) {
    curr_slice = which(meta$slice == slice)
    pixel_barcodes = rownames(meta)[curr_slice]
    pos = Matrix::Matrix(0, nrow = length(curr_slice), ncol = 2)
    pos[,1] = meta[pixel_barcodes,]$x
    pos[,2] = meta[pixel_barcodes,]$y
    colnames(pos) = c('x', 'y')
    rownames(pos) = meta[pixel_barcodes,]$pixel_ID
    pos = as.matrix(pos)
    cellTypeProp = as.matrix(deconvResult[pixel_barcodes,])
    rownames(cellTypeProp) = pixel_barcodes
    title = paste("Slice (", slice, ")", sep = "")
    if (is.na(topic)) {
      a = visAllTopics(cellTypeProp, pos, r=r,topicCols = color, plotTitle = title,
                       showLegend = F,verbose=verbose)
    } else {
      a = STdeconvolve::vizTopic(cellTypeProp, pos, topic=topic, plotTitle = title,
                                 showLegend=F, size = r)
    }
    plot_list <- append(plot_list, list(a))
  }
  return(plot_list)
}



#' Helper function to plot the grid of deconvolution results
#'
#' @description Given the list of plots generated before, now we want to organize
#' them nicely and plot them
#'
#' @param plot_lst Previous plot list
#' @param num_row number of rows for the grid
#' @param num_col number of columns for the grid
#' @param legend_grob object that will become the legend (Default: NA)
#' @param titel title of the graph (Default: "Cell Type Proportion per Pixel")
#'
#' @return Final plot for the deconvolution result
#'
#' @noRd
#'
#' @keywords internal

plot_a_grid <- function(plot_lst, num_row, num_col,
                        legend_grob=NA,
                        title = "Cell Type Proportion per Pixel") {
  plot_grid <- do.call(gridExtra::arrangeGrob,
                       c(plot_lst, ncol = num_col, nrow = num_row,top = title))
  if (length(legend_grob) == 1 && is.na(legend_grob)) {
    combined <- gridExtra::arrangeGrob(plot_grid,
                            ncol = 1,
                            widths = c(1))
  } else {
    combined <- gridExtra::arrangeGrob(plot_grid, legend_grob,
                            ncol = 2,
                            widths = c(4, 1))  # Adjust legend width as needed
  }
  grid::grid.newpage()
  grid::grid.draw(combined)
}






#' Visualize pixel proportions of a single cell-type.
#'
#' @description Note: visualizes all cell-types in theta at once
#' (could be individual cell-types or cell-type-clusters) so for accuracy of
#' the proportions of each cell-type in a pixel, the row (pixel) should sum to 1.
#' This is borrowed from STdeconvolve that is cited below, with some minor
#' adjustments.
#'
#' @param theta document (pixel) x cell-type proportion matrix
#' @param pos position of pixels, as data.frame with x and y columns
#' @param topicOrder order of topics in theta to visualize as a numeric vector
#'     and same length as topicCols (default: seq(ncol(theta)))
#' @param topicCols Vector of colors for each of the cell-types to be visualized.
#'     Same length and order as topicOrder (default: rainbow(ncol(theta)))
#' @param groups Indicates color of the scatterpie strokes (borders) with the goal of coloring them
#'     by their assigned group. This can be a vector or factor indicating the group of each
#'     pixel. Needs to be in the same order as the pixel rows in "theta" (default: NA)
#' @param group_cols Color labels for the groups. Can be a vector or factor. (default: NA)
#' @param r Radius of the scatterpie circles. Adjust based on positions of pixels (default: max(0.4, max(pos)/nrow(pos)*4))
#' @param lwd Width of lines of the pie charts. Increasing helps visualize
#'     group_cols if being used (default: 0.01)
#' @param showLegend Boolean to show the legend indicating cell-types and their color
#' @param plotTitle add title to the resulting plot (default: NA)
#' @param overlay raster image of an H&E tissue (for example) to plot the scatterpies on top of
#'     (default: NA)
#'
#' @return a plot of scatterpies, where each scatterpie represents
#'     a pixel in space based on the x,y coordinates and the components
#'     represent the proportion of each cell-type at that pixel.
#' @source STdeconvolve (https://github.com/JEFworks-Lab/STdeconvolve/blob/package/R/plot.R)
#'
#' @noRd
#'
#' @keywords internal

visAllTopics <- function(theta, pos,
                         topicOrder=seq(ncol(theta)),
                         topicCols=rainbow(ncol(theta)),
                         groups = NA,
                         group_cols = NA,
                         r = max(0.4, max(pos)/nrow(pos)*4),
                         lwd = 0.01,
                         showLegend = TRUE,
                         plotTitle = "Title",
                         overlay = NA,
                         verbose = T) {

  ## check that theta and pos are either data.frames or matrices
  if( !is.matrix(theta) & !is.data.frame(theta) ){
    stop("`theta` must be a matrix or data.frame.")
  }
  if( !is.matrix(pos) & !is.data.frame(pos) ){
    stop("`pos` must be a matrix or data.frame with exactly 2 columns named `x` and `y`.")
  }

  if( (any(!colnames(pos) %in% c("x", "y")) == TRUE) | (dim(pos)[2] != 2) ){
    stop("`pos` must have exactly 2 columns named `x` and `y`.")
  }

  # pixel cell-type distribution reordered based on topicOrder
  theta_ordered <- theta[, topicOrder]
  theta_ordered <- as.data.frame(theta_ordered)
  colnames(theta_ordered) <- paste0("X", colnames(theta_ordered))

  # ensure that `theta` and `pos` pixel rownames maintain same order
  # after the merge so as to not mess up the order of `groups`
  # if provided
  # make sure only using the shared pixels
  pixels <- intersect(rownames(theta_ordered), rownames(pos))
  pixels <- rownames(theta_ordered)[which(rownames(theta_ordered) %in% pixels)]

  # add columns "x", "y" with document positions from `pos`
  theta_ordered_pos <- merge(data.frame(theta_ordered),
                             data.frame(pos), by=0)
  rownames(theta_ordered_pos) <- theta_ordered_pos[,"Row.names"]
  ## make sure pixels in the original order before the merge
  theta_ordered_pos <- theta_ordered_pos[pixels,]

  # first column after merge is "Row.names", last two are "x" and "y"
  # problem is that data frame will replace "-" and " " with "."
  topicColumns <- colnames(theta_ordered_pos)[2:(dim(theta_ordered_pos)[2]-2)]

  # color of piechart groups (lines of piechart):
  if (is.na(groups[1])) {
    groups <- rep("0", dim(theta_ordered_pos)[1])
    theta_ordered_pos$Pixel.Groups <- groups
  } else {
    theta_ordered_pos$Pixel.Groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
    group_cols <- c("0" = "gray")
  }
  if (verbose) {
    message("Plotting scatterpies for ", dim(theta_ordered_pos)[1], " pixels with ", length(topicColumns),
            " cell-types...this could take a while if the dataset is large.", "\n")
  }

  colname = gsub(" ", ".", colnames(theta))

  if (!is.na(overlay[1])){
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = 0:dim(overlay)[2], y = 0:dim(overlay)[1])) +
      ggplot2::coord_equal(xlim = c(0,dim(overlay)[2]), ylim = c(0, dim(overlay)[1]), expand = FALSE) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12, colour = "black"),
        legend.title = ggplot2::element_text(size = 12, colour = "black")
      ) +
      ggplot2::annotation_raster(overlay, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = Pixel.Groups),
                                  lwd = lwd,
                                  data = theta_ordered_pos,
                                  cols = topicColumns,
                                  legend_name = "CellTypes") +
      ggplot2::scale_fill_manual(values = setNames(topicCols, paste0("X", colname))) +
      ggplot2::scale_color_manual(values = group_cols)
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12, colour = "black"),
        legend.title = ggplot2::element_text(size = 12, colour = "black")
      ) +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = Pixel.Groups),
                                  lwd = lwd,
                                  data = theta_ordered_pos,
                                  cols = topicColumns,
                                  legend_name = "CellTypes") +
      ggplot2::scale_fill_manual(values = setNames(topicCols, paste0("X", colname))) +
      ggplot2::scale_color_manual(values = group_cols)
  }

  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }

  p <- p + ggplot2::coord_equal()

  return(p)
}

#' Helper function to retrieve legend object
#'
#' @description Retrieve legend object from ggplot2 plot
#'
#' @param p A ggplot2 plot
#'
#' @return A legend object if it exists
#'
#' @noRd
#'
#' @keywords internal

get_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(legend_index) > 0) {
    return(g$grobs[[legend_index]])
  } else {
    stop("No legend found in the plot.")
  }
}

#' Create legend object for plotting deconvolution result
#'
#' @description Create a pseudo ggplot2 plot using the provided label and
#' color, then retrieving only the legend object and return
#'
#' @param legend_labels Labels for the cell types
#' @param legend_colors Colors for the cell types
#'
#' @return A legend object if it exists
#'
#' @examples
#'
#' data(mouse_hypothalamus)
#' color = c('Astrocyte'='red', 'Endothelial'='orange',
#'           'Ependymal'='black', 'Excitatory'='blue',
#'           'Inhibitory'='green', 'Microglia'='purple',
#'           'OD Immature'='yellow', 'OD Mature'='brown',
#'           'Pericytes'='cyan')
#' legend_labels <- names(color)
#' legend_colors <- as.vector(color)
#' legend_grob <- getLegendGrob(legend_labels, legend_colors)
#' plotDeconvRes(mouse_hypothalamus$ground_truth_deconv,
#'               mouse_hypothalamus$spatial_meta,
#'               color = color,
#'               numRow = 1,
#'               numCol = 2,
#'               legend_grob = legend_grob)
#'
#' @export

getLegendGrob <- function(legend_labels, legend_colors) {
  legend_plot <- ggplot2::ggplot(data.frame(x = 1, y = 1, type = legend_labels),
                                 ggplot2::aes(x = x, y = y, fill = type)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = legend_colors, name = "Cell Types") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")
  legend_grob <- get_legend(legend_plot)
  return(legend_grob)
}

