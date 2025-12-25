#' Assign original cell type identity to estimated cell types
#'
#' @description Assign identity to cell types whose estimated gene expression
#' profile is similar (in terms of correlation) to their original provided gene
#' expression profile.
#'
#' @param estGeneExp unscaled estimated cell type gene expression profile
#' @param provGeneExp provided cell type gene expression profile
#' @param deconvRes unscaled deconvolution result
#' @param unknownIdentity cell types with unknown identities, for example placeholders
#'  (Default: character(0))
#' @param corMethod correlation method (Default: "pearson")
#' @param threshold correlation threshold for assigning cell type identity (Default: 0.8)
#' @param verbose Default: T
#'
#' @return A list containing deconvolution result and gene expression profile for
#' cell types with identity assigned, as well as cell types with unknown identity
#'
#' @examples
#'
#'
#' @export


identityAssignment <- function(estGeneExp, provGeneExp, deconvRes,
                               unknownIdentity = character(0),
                               corMethod = "pearson",
                               threshold = 0.8,
                               verbose = T) {

  # step 1 of interpretation
  est_gene_exp = estGeneExp/rowSums(as.matrix(estGeneExp))
  prov_gene_exp = provGeneExp/rowSums(as.matrix(provGeneExp))
  deconv_res = deconvRes/rowSums(deconvRes)
  est_ct_corr = cor(t(as.matrix(est_gene_exp)),
                    t(as.matrix(prov_gene_exp)),
                    method = corMethod)
  known_identity = c()
  unknown_identity = unknownIdentity


  # first figure out how many cell types have identity assigned
  for (ct in colnames(est_ct_corr)) {
    if (ct %in% unknown_identity) {
      next
    }
    if (est_ct_corr[ct, ct] > threshold) {
      if (verbose) {
        message(paste(ct, "is assigned to itself"))
      }
      known_identity = c(known_identity, ct)
    } else {
      unknown_identity = c(unknown_identity, ct)
    }
  }
  return(list(assignedDeconvRes = deconvRes[,known_identity],
              assignedGeneExp = estGeneExp[known_identity,],
              unknownIdentity = unknown_identity))
}



#' Assign non-original cell type identity to estimated cell types
#'
#' @description Assign identity to cell types whose estimated gene expression
#' profile is similar (in terms of correlation) to provided gene expression profile
#' of any other cell types
#'
#' @param estGeneExp unscaled estimated cell type gene expression profile
#' @param provGeneExp provided cell type gene expression profile
#' @param deconvRes unscaled deconvolution result
#' @param unknownIdentity cell types with unknown identities so far
#' @param assignedDeconvRes unscaled estimated deconvolution result of cell types
#' with identity already assigned
#' @param assignedGeneExp unscaled estimated gene expression profile of cell types
#' with identity already assigned
#' @param corMethod correlation method (Default: "pearson")
#' @param threshold correlation threshold for assigning cell type identity (Default: 0.8)
#' @param verbose Default: T
#'
#' @return A list containing deconvolution result and gene expression profile for
#' cell types with identity assigned, as well as cell types with unknown identity
#'
#' @examples
#'
#'
#' @export

identityAssignmentOther <- function(estGeneExp, provGeneExp, deconvRes,
                                    unknownIdentity, assignedDeconvRes,
                                    assignedGeneExp,
                                    corMethod = "pearson",
                                    threshold = 0.8,
                                    verbose = T) {
  # step 2 of interpretation, we want to match to other cell types
  est_gene_exp = estGeneExp/rowSums(as.matrix(estGeneExp))
  prov_gene_exp = provGeneExp/rowSums(as.matrix(provGeneExp))
  deconv_res = deconvRes/rowSums(deconvRes)
  # index goes: est_ct_corr[estimated, provided]
  est_ct_corr = cor(t(as.matrix(est_gene_exp)),
                    t(as.matrix(prov_gene_exp)),
                    method = corMethod)
  unknown_identity = unknownIdentity
  known_identity = setdiff(rownames(estGeneExp), unknown_identity)

  assigned_deconv_res = assignedDeconvRes
  assigned_gene_exp = assignedGeneExp
  labelled_ct = colnames(assignedDeconvRes)

  for (ct in unknown_identity) {
    cor_vec = est_ct_corr[ct,]
    cor_vec = cor_vec[names(cor_vec) != ct]

    i <- which.max(cor_vec)        # index of maximum correlation
    max_corr <- cor_vec[i]
    max_ct  <- names(cor_vec)[i]
    if (max_corr > threshold) {
      known_identity = c(known_identity, ct)
      if (verbose) {
        message(paste(ct, "is assigned to", max_ct))
      }

      if (max_ct %in% labelled_ct) {
        assigned_deconv_res[,max_ct] = assigned_deconv_res[,max_ct] + deconvRes[,ct]
        assigned_gene_exp[max_ct,] = assigned_gene_exp[max_ct,] + estGeneExp[ct,]
      } else {
        assigned_deconv_res = cbind(assigned_deconv_res, deconvRes[,ct])
        colnames(assigned_deconv_res)[ncol(assigned_deconv_res)] = max_ct

        assigned_gene_exp = rbind(assigned_gene_exp, estGeneExp[ct,])
        rownames(assigned_gene_exp)[nrow(assigned_gene_exp)] = max_ct
      }
    }
  }
  return(list(assignedDeconvRes = assigned_deconv_res,
              assignedGeneExp = estGeneExp[known_identity,],
              unknownIdentity = setdiff(rownames(estGeneExp), known_identity)))

}

#' For the final unknown cell types, perform clustering
#'
#' @description For the final unknown cell types, perform clustering and report
#' final interpreted deconvolution result and gene expression profile
#'
#' @param estGeneExp unscaled estimated cell type gene expression profile
#' @param deconvRes unscaled deconvolution result
#' @param unknownIdentity cell types with unknown identities so far
#' @param assignedDeconvRes unscaled estimated deconvolution result of cell types
#' with identity already assigned
#' @param assignedGeneExp unscaled estimated gene expression profile of cell types
#' with identity already assigned
#' @param corMethod correlation method (Default: "pearson")
#' @param threshold correlation threshold for assigning cell type identity (Default: 0.7)
#' @param verbose Default: T
#'
#' @return A list containing deconvolution result and gene expression profile for
#' cell types with identity assigned, as well as cell types with unknown identity
#'
#' @examples
#'
#'
#' @export

unknownCluster <- function(estGeneExp, deconvRes,
                           unknownIdentity, assignedDeconvRes,
                           assignedGeneExp,
                           corMethod = "pearson",
                           threshold = 0.7,
                           verbose = T) {
  # step 2 of interpretation, we want to match to other cell types
  est_gene_exp = estGeneExp/rowSums(as.matrix(estGeneExp))
  deconv_res = deconvRes/rowSums(as.matrix(deconvRes))
  # index goes: est_ct_corr[estimated, provided]
  est_ct_corr = cor(t(as.matrix(est_gene_exp)),
                    t(as.matrix(est_gene_exp)),
                    method = corMethod)
  unknown_identity = unknownIdentity
  known_identity = setdiff(rownames(estGeneExp), unknown_identity)

  assigned_deconv_res = assignedDeconvRes
  assigned_gene_exp = assignedGeneExp
  unassigned_deconv_res = deconvRes[,unknown_identity]
  unassigned_gene_exp = estGeneExp[unknown_identity,]

  #labelled_ct = colnames(assignedDeconvRes)
  abundance = sort(colSums(deconv_res[,unknown_identity]), decreasing=F)

  clustered_ct = c()
  final_unknowns = c()

  for (ct in names(abundance)) {
    clustered_ct = c(clustered_ct, ct)
    remaining_ct = setdiff(names(abundance), clustered_ct)
    if (length(remaining_ct) == 0) {
      final_unknowns = c(final_unknowns, ct)
      break
    }
    cor_vec = est_ct_corr[ct, remaining_ct]

    i <- which.max(cor_vec)        # index of maximum correlation
    max_corr <- cor_vec[i]
    max_ct  <- names(cor_vec)[i]

    if (max_corr > threshold) {
      # then we merge
      if (verbose) {
        message(paste("Estimated", ct, "is merged with estimated", max_ct))
      }
      unassigned_deconv_res[,max_ct] = unassigned_deconv_res[,max_ct] + deconvRes[,ct]
      unassigned_gene_exp[max_ct,] = unassigned_gene_exp[max_ct,] + estGeneExp[ct,]
    } else {
      final_unknowns = c(final_unknowns, ct)
    }
  }

  unassigned_deconv_res = unassigned_deconv_res[,final_unknowns]
  unassigned_gene_exp = unassigned_gene_exp[final_unknowns,]
  for (i in 1:length(final_unknowns)) {
    if (verbose) {
      message(paste("Estimated",final_unknowns[i],"becomes Unknown",i))
    }
  }
  colnames(unassigned_deconv_res) = paste("Unknown",1:length(final_unknowns))
  rownames(unassigned_gene_exp) = paste("Unknown",1:length(final_unknowns))

  finalDeconvRes = cbind(assignedDeconvRes, unassigned_deconv_res)
  finalGeneExp = rbind(assignedGeneExp, unassigned_gene_exp)

  finalDeconvRes <- finalDeconvRes[, order(colnames(finalDeconvRes))]
  finalGeneExp <- finalGeneExp[order(rownames(finalGeneExp)),]

  return(list(finalDeconvRes = finalDeconvRes,
              finalGeneExp = finalGeneExp))
}











