library(ggplot2)
library(STdeconvolve)
library(Matrix)
library(gridExtra)
library(pheatmap)
library(bayess)

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
        legend.title = ggplot2::element_text(size = 12, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5)
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
        legend.title = ggplot2::element_text(size = 12, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5)
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

get_legend <- function(p) {
  g <- ggplotGrob(p)
  legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(legend_index) > 0) {
    return(g$grobs[[legend_index]])
  } else {
    stop("No legend found in the plot.")
  }
}

plot_a_grid <- function(plot_grid, legend_grob=NA, no_legend=FALSE, output_file=NULL) {
  if (no_legend) {
    combined <- gridExtra::arrangeGrob(plot_grid,
                                       ncol = 1,
                                       widths = c(3))
  } else {
    combined <- gridExtra::arrangeGrob(plot_grid, legend_grob,
                                       ncol = 2,
                                       widths = c(4, 1))  # Adjust legend width as needed
  }

  # If output_file is provided, save as SVG
  if (!is.null(output_file)) {
    # Open SVG device
    svg(filename = output_file, width = 8, height = 6)  # adjust width/height as needed
    grid.draw(combined)
    dev.off()  # close the device
  } else {
    # Just draw in current device
    grid.newpage()
    grid.draw(combined)
  }
}

reorder_celltypes <- function(reordering, deconvolvedProp, deconvolvedGexp) {
  matched.deconProp = Matrix(0, nrow = dim(deconvolvedProp)[1],
                             ncol = dim(deconvolvedProp)[2])
  matched.deconGexp = Matrix(0, nrow = dim(deconvolvedGexp)[1],
                             ncol = dim(deconvolvedGexp)[2])
  for (i in 1:(length(reordering1))) {
    matched.deconProp[,i] = deconvolvedProp[,reordering[i]]
    matched.deconGexp[i,] = deconvolvedGexp[reordering[i],]
  }
  matched.deconProp = as.matrix(matched.deconProp)
  matched.deconGexp = as.matrix(matched.deconGexp)

  rownames(matched.deconProp) = rownames(deconvolvedProp)
  rownames(matched.deconGexp) = rownames(deconvolvedGexp)
  colnames(matched.deconGexp) = colnames(deconvolvedGexp)
  return(list(reorder_deconProp = matched.deconProp,
              reorder_deconGexp = matched.deconGexp))

}






annotateCellTypesGSEA <- function(beta, gset, qval=0.05, ...) {
  # Modified based on STdeconvolve
  results <- list()
  top.pos.enrich <- c()

  for (i in seq(nrow(beta))){
    celltype <- i
    vals <- sort(beta[celltype,], decreasing=TRUE)

    gsea.results <- fgsea::fgsea(stats=vals, pathways=gset, scoreType = "pos", ...)

    # filter for top hits
    gsea.sig <- gsea.results[gsea.results$padj < qval,]

    if(nrow(gsea.sig) > 1) {
      gsea.sig <- gsea.sig[order(gsea.sig$padj, -sapply(gsea.sig$leadingEdge, length), -gsea.sig$ES), ]
    }


    results[[ rownames(beta)[celltype] ]] <- gsea.sig

    ## the top entry that is also positiviely enriched in the txn profile is predicted to be the best matching
    gsea.sig.pos <- gsea.sig[which(gsea.sig$ES > 0), ]
    top.pos.enrich <- append(top.pos.enrich, gsea.sig.pos[1]$pathway)

  }

  #names(top.pos.enrich) <- rownames(beta)

  return(list(results=results,
              predictions=top.pos.enrich))

}






getCellTypeMarkers <- function(gt_reference, top_genes = 20) {
  cellTypeMarkers <- list()
  gexp <- as.matrix(gt_reference*100)

  for (i in seq(length(rownames(gexp)))){
    celltype <- i
    ## log2FC relative to other cell-types
    ## highly expressed in cell-type of interest
    highgexp <- names(sort(gexp[celltype,], decreasing = TRUE))[1:top_genes]
    ## high log2(fold-change) compared to other deconvolved cell-types and limit
    ## to the top 200
    log2fc <- sort(
      log2(gexp[celltype,highgexp]/colMeans(gexp[-celltype,highgexp])),
      decreasing=TRUE)

    ## for gene set of the ground truth cell-type, get the genes
    ## with log2FC > 1 (so FC > 2 over the mean exp of the other cell-types)
    markers <- names(log2fc[log2fc > 1])
    cellTypeMarkers[[ rownames(gexp)[celltype] ]] <- markers
  }
  return(cellTypeMarkers)
}




cor_matrix_plot <- function(gt_ref, deconGexp, title = "Transcriptional correlation",
                            type = "pearson", rowlab = "Deconvolved cell-types",
                            collab = "Ground truth cell-types") {
  cormatrix = cor(t(gt_ref), t(deconGexp), method = type)

  #pairs <- lsatPairs(cormatrix)
  #cormatrix2 = cor(t(gt_ref), t(deconGexp)[,pairs$colsix], method=type)
  #colnames(cormatrix2) <- paste0("X", colnames(cormatrix2))
  p = correlationPlot(mat = cormatrix,
                      rowLabs = rowlab, # aka x-axis, and rows of matrix
                      colLabs = collab, # aka y-axis, and columns of matrix
                      title = title, annotation = TRUE) +

    ## this function returns a `ggplot2` object, so can add additional aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))
  p$layers[[2]]$aes_params$size <- 2.5
  return(p)
}

cor_matrix_plot_gt <- function(gt_ref, title = "Transcriptional correlation",
                               type = "pearson", axis_label = "Ground truth cell-types") {
  cormatrix = cor(t(gt_ref), t(gt_ref), method=type)
  pairs <- lsatPairs(cormatrix)
  cormatrix2 = cor(t(gt_ref), t(gt_ref)[,pairs$colsix], method=type)
  p = correlationPlot(mat = cormatrix2,
                      colLabs = axis_label,
                      rowLabs = axis_label,
                      title = title,
                      annotation = TRUE) +

    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0, size = 15),  # smaller x-axis labels
      axis.text.y = ggplot2::element_text(size = 15),                        # smaller y-axis labels
      axis.title.x = ggplot2::element_text(size = 15),                       # smaller x-axis title
      axis.title.y = ggplot2::element_text(size = 15),                       # smaller y-axis title
      plot.title = ggplot2::element_text(size = 15)                         # smaller main title
    )

  p$layers[[2]]$aes_params$size <- 2.75

  return(p)
}



is_number <- function(x) {
  !is.na(as.numeric(x))
}



plot_pixels <- function(gt, meta, verbose=T, topic=NA, color_prov=NA) {
  normalized_gt_molecule = (gt/rowSums(gt))
  if (!is.na(color_prov)) {
    print("Provided color")
    cols = color_prov
  } else {
    cols = c('Astrocyte'='red', 'Endothelial'='orange',
             'Ependymal'='black', 'Excitatory'='blue',
             'Inhibitory'='green', 'Microglia'='purple',
             'OD Immature'='yellow', 'OD Mature'='brown',
             'Pericytes'='cyan')
  }

  plots_list <- list()
  for (each_bregma in unique(meta$bregma)) {
    curr_bregma = which(meta$bregma == each_bregma)
    center_x = 0.5*(meta$bottom_left_x[curr_bregma] + meta$top_right_x[curr_bregma])
    center_y = 0.5*(meta$bottom_left_y[curr_bregma] + meta$top_right_y[curr_bregma])
    pos = Matrix(0, nrow = length(curr_bregma), ncol = 2)
    colnames(pos) = c('x', 'y')
    pos[,1] = center_x
    pos[,2] = center_y
    rownames(pos) = meta$pixel_ID[curr_bregma]
    pos = as.matrix(pos)
    title = paste("Bregma (",each_bregma,")",
                  sep="")
    if (is.na(topic)) {
      a = visAllTopics(as.matrix(normalized_gt_molecule[curr_bregma,]), pos, r=40,
                       topicOrder = seq(length(cols)), topicCols = cols,
                       showLegend =F,plotTitle = title,verbose=verbose)
    } else {
      a = vizTopic(as.matrix(normalized_gt_molecule[curr_bregma,]), pos,
                   topic=topic, showLegend = F, plotTitle=title)
    }

    #plot(a)
    plots_list <- append(plots_list, list(a))
  }
  return(plots_list)
}


plot_est_pixels <- function(deconv_result, meta, cell_type_names, verbose=T,
                            new_color = F, topic=NA, new_cols=NA) {
  cols = c('Astrocyte'='red', 'Endothelial'='orange',
           'Ependymal'='black', 'Excitatory'='blue',
           'Inhibitory'='green', 'Microglia'='purple',
           'OD Immature'='yellow', 'OD Mature'='brown',
           'Pericytes'='cyan')
  if (new_color) {
    cols = new_cols
  }

  plot_list = list()
  for (each_bregma in unique(meta$bregma)) {
    curr_bregma = which(meta$bregma == each_bregma)
    center_x = 0.5*(meta$bottom_left_x[curr_bregma] + meta$top_right_x[curr_bregma])
    center_y = 0.5*(meta$bottom_left_y[curr_bregma] + meta$top_right_y[curr_bregma])
    pos = Matrix(0, nrow = length(curr_bregma), ncol = 2)
    pos[,1] = center_x
    pos[,2] = center_y
    colnames(pos) = c('x', 'y')
    rownames(pos) = meta$pixel_ID[curr_bregma]
    pos = as.matrix(pos)

    cellTypeProp = (deconv_result$gamma/rowSums(deconv_result$gamma))[curr_bregma,]
    rownames(cellTypeProp) = meta$pixel_ID[curr_bregma]
    colnames(cellTypeProp) = cell_type_names
    title = paste("Bregma (", each_bregma, ")", sep = "")

    if (is.na(topic)) {
      a = visAllTopics(cellTypeProp, pos, r=40,topicCols = cols, plotTitle = title,
                       showLegend = F,verbose=verbose)
    } else {
      a = vizTopic(cellTypeProp, pos, topic=topic, plotTitle = title, showLegend=F,
                   verbose=verbose)
    }

    plot_list <- append(plot_list, list(a))
  }
  return(plot_list)
}

create_zero_col <- function(mat, specific_name) {
  # Convert sparse matrix to a regular dense matrix
  mat <- as.matrix(mat)

  # Convert to a data frame to allow character column indexing
  mat <- as.data.frame(mat)

  # Add a new column of 0s with the correct name
  mat[[specific_name]] <- 0

  # Sort columns alphabetically
  mat <- mat[, order(colnames(mat))]

  return(as.matrix(mat))
}




find_outlier <- function(data_vector) {
  # This function finds non-outliers of a data vector by 1.5 IQR above Q3 or
  # below Q1 respectively
  Q1 <- quantile(data_vector, 0.25)
  Q3 <- quantile(data_vector, 0.75)
  IQR_value <- Q3 - Q1

  # Compute lower and upper bounds
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value

  # Identify indices of non-outliers
  non_outlier_indices <- which(data_vector >= lower_bound & data_vector <= upper_bound)
  return(non_outlier_indices)
}


remove_cellType_RCTD <- function(ct1, ct2, celltypeList, corpus_count, umi,
                                 puck, output_name) {
  miss_ct1_idx = which(celltypeList != ct1)
  miss_ct2_idx = which(celltypeList != ct2)
  miss_ct_idx = intersect(miss_ct1_idx, miss_ct2_idx)
  count = corpus_count[,miss_ct_idx]
  ct_remain = celltypeList[miss_ct_idx]
  ct_remain <- factor(ct_remain, levels = unique(ct_remain))
  UMI_remain = umi[miss_ct_idx]
  curr_ref = Reference(count, ct_remain, UMI_remain)

  myRCTD <- create.RCTD(puck, curr_ref, max_cores = 1, UMI_min = 1, counts_MIN = 1)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD@results$weights, output_name)
}


remove_ref_cellTypes <- function(ct_list, reference, pw, pw_output,new_ref_output,
                                 randomize=F, randomize_param=50) {
  # reference should be cell type x gene

  # firstly, let's check if ct_list are within reference
  for (ct in ct_list) {
    if (!(ct %in% rownames(reference))) {
      print(paste("Cell Type", ct, "not a valid cell type"))
      return(NULL)
    }

  }

  num_celltype = dim(reference)[1]
  fully_trust_reference = rep(1, num_celltype)
  for (i in 1:num_celltype) {
    if (rownames(reference)[i] %in% ct_list) {
      # it is one of the cell types that we ignore
      fully_trust_reference[i] = pw
    }
  }

  ref_copy = reference
  if (randomize) {
    ref_copy[ct_list,] = rdirichlet(n=length(ct_list),
                      par = rep(randomize_param,dim(reference)[2]))
  } else {
    ref_copy[ct_list,] = 1/dim(reference)[2]
  }

  saveRDS(ref_copy, new_ref_output)
  saveRDS(fully_trust_reference, pw_output)
  return("Success")
}


visualize_reference <- function(ref, want_log=F, cluster_col=T, cluster_row=F,
                                title = "Gene Expression Heatmap") {


  if (want_log) {
    new_ref = log(ref)
  } else {
    new_ref = ref
  }

  min_val <- min(new_ref)
  max_val <- max(new_ref)

  heatmap_plot <- pheatmap(new_ref,
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





compute_RMSE <- function(mat1, mat2, col_check = F, renorm = F) {

  if (col_check) {
    common_cols = sort(intersect(colnames(mat1), colnames(mat2)))
    mat1 = mat1[,common_cols]
    mat2 = mat2[,common_cols]
  }

  if (renorm) {
    mat1 = mat1/rowSums(mat1)
    mat2 = mat2/rowSums(mat2)
    mat1[is.nan(mat1)] = 1/length(colnames(mat1))
    mat2[is.nan(mat2)] = 1/length(colnames(mat2))

  }

  diff = rowSums((mat1 - mat2)^2)
  sqrt(diff/dim(mat1)[2])
}

KL_divergence <- function(ground_truth, estimate) {
  # Calculates KL divergence for categorical distribution
  # Need to ensure estimate has no 0s
  log_p_div_q = log(ground_truth/estimate)
  log_p_div_q[is.infinite(log_p_div_q)] = 0
  p_log_p_div_q = ground_truth * log_p_div_q
  return(rowSums(p_log_p_div_q))
}



violin_plot <- function(df, x_col, y_col, remove_outlier = F,
                        title = "RMSE comparison between models") {
  if (remove_outlier) {
    # Calculate the IQR for the y-variable (outlier detection)
    Q1 <- quantile(df[[y_col]], 0.25, na.rm = TRUE)
    Q3 <- quantile(df[[y_col]], 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1

    # Remove outliers using the 1.5 * IQR rule
    df_filtered <- df[df[[y_col]] >= (Q1 - 1.5 * IQR) & df[[y_col]] <= (Q3 + 1.5 * IQR), ]
  } else {
    df_filtered = df
  }


  # Create the violin plot with filtered data
  p1 = ggplot(df_filtered, aes_string(x = x_col, y = y_col, fill = x_col)) +
    geom_violin(trim = FALSE) +
    labs(title = title,
         x = "Model", y = y_col) +
    stat_summary(fun = "mean", geom = "point", size = 3, color = "black", shape = 18, show.legend = FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p1)
}



reorder_vector <- function(vec, values_to_front) {
  # Reorder a vector vec such that front_values appear in the front
  front_values <- vec[vec %in% values_to_front]
  remaining_values <- vec[!(vec %in% values_to_front)]
  # Combine them
  reordered_vec <- c(front_values, remaining_values)
  return(reordered_vec)
}


vis_ref_for_miss_celltype <- function(deconv_result, ref_mat, ct, orig_order,
                                      celltypes_front = c(), want_log = F,
                                      cluster_col = T) {
  gene_exp = deconv_result
  rownames(gene_exp) = paste("Estimated", rownames(gene_exp))

  if (length(celltypes_front) == 0) {
    # no predefined order
    gene_exp_together = rbind(ref_mat, gene_exp)[reorder_vector(orig_order,
        c(ct, paste("Estimated", ct), "Inhibitory", "Estimated Inhibitory")),]
  } else {
    gene_exp_together = rbind(ref_mat, gene_exp)[reorder_vector(orig_order,
                                                          celltypes_front),]
  }


  gene_exp_plot = visualize_reference(gene_exp_together, want_log = want_log,
                                      cluster_col = cluster_col,
                        title = paste("Gene Expression after removing", ct))
  return(gene_exp_plot)
}




dirichlet_var <- function(mat) {
  # matrix, where each row represents parameters of a dirichlet distribution
  alpha_0 = rowSums(mat)
  alpha_tilde = mat/alpha_0
  return(alpha_tilde*(1-alpha_tilde)/(alpha_0 + 1))
}
