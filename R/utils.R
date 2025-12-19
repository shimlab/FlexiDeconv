
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
#' data(mouse_hypothalamus)
#' new_reference <- appendPlaceholder(mouse_hypothalamus$reference, numPlaceholder=1)
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












