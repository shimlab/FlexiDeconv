
#' Append placeholder cell types onto the reference
#'
#' @description Given current cell type x gene reference matrix, append placeholder
#' cell types with number specified by the user, default to be 1. Note that placeholder
#' will be generated using dirichlet distribution with uniform mean
#'
#' @param reference original reference matrix
#' @param numPlaceholder number of placeholder cell types wanted (Default: 1)
#'
#' @return New reference matrix with dimension (cell type + placeholder) x gene
#'
#' @examples
#'
#' data(mouseHypothalamus)
#' new_reference <- appendPlaceholder(mouseHypothalamus$reference, numPlaceholder=1)
#'
#' @export



appendPlaceholder <- function(reference, numPlaceholder = 1) {
  if (numPlaceholder == 0) {
    return(reference)
  }
  num_gene = ncol(reference)
  new_reference = rbind(reference,
                        bayess::rdirichlet(n=numPlaceholder,
                                           par = rep(5*num_gene,num_gene)))
  rownames(new_reference) = c(rownames(reference),
                              paste("Placeholder",1:numPlaceholder))
  return(new_reference)
}




#' For provided spatial data and scRNA-seq data, find intersecting gene set
#'
#' @description For provided spatial data and scRNA-seq data, find intersecting gene set
#'
#' @param spatial matrix representing spatial data with dimension pixel x gene
#' @param scRNAseq matrix representing scRNAseq data with dimension cell x gene
#'
#' @return spatial and scRNA-seq data after selecting only common gene set
#'
#' @examples
#'
#' data(scRNAseq_to_reference)
#' common_geneset <- commonGeneFilter(scRNAseq_to_reference$spatial, scRNAseq_to_reference$scRNAseq)
#'
#' @export



commonGeneFilter <- function(spatial, scRNAseq) {

  common_geneset = intersect(colnames(spatial), colnames(scRNAseq))
  return(list(spatial = spatial[,common_geneset],
              scRNAseq = scRNAseq[,common_geneset]))
}



#' Remove redundant genes
#'
#' @description For given spatial data and scRNAseq data, remove genes that have
#' 0 counts across all pixels and all cells.
#'
#' @param spatial matrix representing spatial data with dimension pixel x gene
#' @param scRNAseq matrix representing scRNAseq data with dimension cell x gene
#'
#' @return spatial and scRNA-seq data after removing redundant genes
#'
#' @examples
#'
#' data(scRNAseq_to_reference)
#' common_geneset <- commonGeneFilter(scRNAseq_to_reference$spatial, scRNAseq_to_reference$scRNAseq)
#' redundant_removal <- redundantGeneFilter(common_geneset$spatial, common_geneset$scRNAseq)
#'
#' @export



redundantGeneFilter <- function(spatial, scRNAseq) {

  totalGeneCountSpatial = colSums(as.matrix(spatial))
  totalGeneCountscRNA = colSums(as.matrix(scRNAseq))

  genes_to_keep = intersect(names(totalGeneCountscRNA[totalGeneCountscRNA>0]),
                            names(totalGeneCountSpatial[totalGeneCountSpatial>0]))
  return(list(spatial = spatial[,genes_to_keep],
              scRNAseq = scRNAseq[,genes_to_keep]))

}


#' Outlier cell removal in scRNAseq data
#'
#' @description For provided scRNAseq data, remove outlier cells. Outlier is
#' defined as:
#' 1) cells with total gene count that is 1.5 IQR above Q3 of all
#' cells within the current cell type.
#' 2) cells with rank x (default: 2) having at least y (default: 5) count.
#'
#' @param scRNAseq matrix representing scRNAseq data with dimension cell x gene
#'
#' @return scRNAseq data after removing outlier cells
#'
#' @examples
#'
#' data(scRNAseq_to_reference)
#' common_geneset <- commonGeneFilter(scRNAseq_to_reference$spatial, scRNAseq_to_reference$scRNAseq)
#'
#' @export



outlierCellRemoval <- function(scRNAseq, celltypes, lowerThresholdRank = 2,
                               lowerThresholdCount = 3, plot = T, plot_density = 5) {

  totalGeneCount = rowSums(scRNAseq)
  cells_left = c()

  celltype_forplot = c()
  totalGeneCount_forplot = c()
  total_celltype_covered = 1
  curr_celltype_covered = 1

  colors = c('red', 'green', 'orange', 'purple', 'cyan')

  for (celltype in unique(celltypes)) {

    curr_cells = names(celltypes)[celltypes == celltype]
    Q3 = quantile(totalGeneCount[curr_cells], 0.75, na.rm = T)
    IQR = IQR(totalGeneCount[curr_cells], na.rm = T)
    upperthreshold = 1.5 * IQR + Q3
    vals <- totalGeneCount[curr_cells]
    within_upper = names(vals)[vals <= upperthreshold]

    sortTotalGeneCount = sort(vals, decreasing = F)
    if (length(vals) < lowerThresholdRank) {
      print(paste("Not enough cells for celltype", celltype))
      next
    }
    rankN <- apply(scRNAseq, 1, function(x) sort(x, decreasing = TRUE)[lowerThresholdRank])
    within_lower = names(rankN[rankN>=lowerThresholdCount])
    cells_satisfy = intersect(within_upper, within_lower)
    cells_left = c(cells_left, cells_satisfy)

    #if (plot == F) {next}
    celltype_forplot = c(celltype_forplot, rep(celltype, length(cells_satisfy)))
    totalGeneCount_forplot = c(totalGeneCount_forplot, totalGeneCount[cells_satisfy])
    if ((total_celltype_covered == length(unique(celltypes))) ||
        (curr_celltype_covered == plot_density)) {
      df = data.frame(cell_type = celltype_forplot,
                      gene_count = totalGeneCount_forplot)

      p = ggplot(df, aes(x = cell_type, y = gene_count, fill = cell_type, color = cell_type)) +
        stat_halfeye(adjust = 0.7, width = 0.6, justification = -0.2, .width = 0) +
        geom_jitter(width = 0.12, alpha = 0.5, size = 1) +
        scale_fill_manual(values = colors[1:curr_celltype_covered]) +
        scale_color_manual(values = colors[1:curr_celltype_covered]) +
        theme_minimal() +
        labs(x = NULL, y = "Total gene count")
      plot(p)

      total_celltype_covered = total_celltype_covered + 1
      curr_celltype_covered = 1
      celltype_forplot = c()
      totalGeneCount_forplot = c()
    } else {
      curr_celltype_covered = curr_celltype_covered + 1
      total_celltype_covered = total_celltype_covered + 1
    }


  }

  return(list(scRNAseq = scRNAseq[cells_left,],
              celltypes = celltypes[cells_left]))
}


#' Form reference based on processed scRNAseq data
#'
#' @description For given scRNAseq data and the cell type annotation for every
#' cell, form cell type x gene refernece.
#'
#' @param scRNAseq matrix representing scRNAseq data with dimension cell x gene
#' @param celltype named vector representing cell type annotation for every cell
#' @param minCell minimum number of cells required to keep current cell type
#' @param emptyFill we replace genes that are not expressed at all in a cell type
#' by a pseudocount (default: 0.1)
#'
#' @return reference with dimension cell type x gene
#'
#' @examples
#'
#' data(scRNAseq_to_reference)
#' common_geneset <- commonGeneFilter(scRNAseq_to_reference$spatial, scRNAseq_to_reference$scRNAseq)
#' redundant_removal <- redundantGeneFilter(common_geneset$spatial, common_geneset$scRNAseq)
#'
#' @export



formReference <- function(scRNAseq, celltype, minCell = 25, emptyFill = 0.1) {

  cell_count = table(celltype)
  celltypes_keep = names(cell_count[cell_count >= minCell])

  ref = Matrix::Matrix(0, nrow = length(celltypes_keep),
                       ncol = ncol(scRNAseq))
  rownames(ref) = sort(celltypes_keep)
  colnames(ref) = colnames(scRNAseq)
  for (curr_celltype in sort(unique(celltype))) {
    if (curr_celltype %in% celltypes_keep) {
      cells = names(celltype[celltype == curr_celltype])
      ref[curr_celltype,] = colSums(scRNAseq[cells,])
    }
  }
  ref[ref == 0] = emptyFill

  return(ref/rowSums(as.matrix(ref)))

}


#' Replace 0s in the reference
#'
#' @description For user-provided reference, replace 0 by 0.1 (default) times the
#' minimum non-zero value in the reference
#'
#' @param reference user-provided reference with dimension cell type x gene
#' @param zeroFillFactor user-configurable factor for determining the factor to
#' multiply the minimum non-zero value by, to replace 0s
#'
#' @return reference with dimension cell type x gene
#'
#' @examples
#'
#' data(reference_ready)
#' ref_replace_zero <- replaceZero(reference_ready, zeroFillFactor = 0.1)
#'
#' @export



replaceZero <- function(reference, zeroFillFactor = 0.1) {

  min_nonzero = min(reference[reference != 0])
  reference[reference == 0] = zeroFillFactor * min_nonzero
  return(reference/rowSums(as.matrix(reference)))

}
