

#' Helper function to calculate the ELBO
#'
#' @description For given reference, gamma, tau, alpha, phi parameter values, as well
#' as cell type counts, pixel counts, gene counts and the spatial transcriptomics
#' data, calculate the ELBO
#'
#' @param ref.p reference matrix, with dimension cell type x gene
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param tau.p Variational Inference parameter tau matrix,
#'      with dimension cell type x gene
#' @param alpha.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param phi.p Variational Inference parameter 3-dimensional phi array,
#'      with dimension pixel x gene x cell type
#' @param ct.count total cell type count
#' @param pixel.count total pixel count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, with dimension pixel x gene
#'
#' @return ELBO value for current parameter set
#'
#' @noRd
#'
#' @keywords internal

ELBO <- function(ref.p, gamma.p, tau.p, alpha.p, phi.p, ct.count,
                 pixel.count, gene.count, count.data) {
  start_time <- proc.time()
  term1 <- t1(ref.p, ct.count, tau.p)
  term2 <- t2(alpha.p, gamma.p, pixel.count)
  term3 <- t3(phi.p, gamma.p, pixel.count, ct.count, count.data, gene.count)
  term4 <- t4(gene.count, pixel.count, ct.count, phi.p, tau.p, count.data)
  term5 <- t5(pixel.count, gamma.p)
  term6 <- t6(pixel.count, ct.count, gene.count, phi.p, count.data)
  term7 <- t7(ct.count, tau.p)
  gc()
  #print(paste("ELBO took:", (proc.time() - start_time)[3], "seconds"))
  term1 + term2 + term3 + term4 - term5 - term6 - term7
}

#' Helper function to calculate first term of the ELBO: E[p(beta)])
#'
#' @description For given reference, tau, cell type counts, calculate first
#' term of the ELBO: E[p(beta)])
#'
#' @param ref.p reference matrix, with dimension cell type x gene
#' @param tau.p Variational Inference parameter tau matrix,
#'      with dimension cell type x gene
#' @param ct.count total cell type count
#'
#' @return Calculated first term E[p(beta)]
#'
#' @noRd

t1 <- function(ref.p, ct.count, tau.p) {
  # E[p(beta)]
  val <- 0
  for (i in 1:ct.count) {
    val <- val + lgamma(sum(ref.p[i,])) - sum(lgamma(ref.p[i,]))
    tau.row <- tau.p[i,]
    tau.vec <- digamma(tau.row) - digamma(sum(tau.row))
    val <- val + (ref.p[i,] - 1) %*% tau.vec
  }
  val
}

#' Helper function to calculate second term of the ELBO: E[log p(theta|alpha)]
#'
#' @description For given gamma, alpha, total pixel counts, calculate second
#' term of the ELBO: E[log p(theta|alpha)]
#'
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param alpha.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param pixel.count total pixel count
#'
#' @return Calculated second term E[log p(theta|alpha)]
#'
#' @noRd
#'
#' @keywords internal

t2 <- function(alpha.p, gamma.p, pixel.count) {
  # E[log p(theta|alpha)]
  val <- pixel.count * (lgamma(sum(alpha.p)) - sum(lgamma(alpha.p)))

  for (i in 1:pixel.count) {
    gamma.row <- gamma.p[i,]
    gamma.vec <- digamma(gamma.row) - digamma(sum(gamma.row))
    val <- val + (alpha.p - 1) %*% gamma.vec
  }
  val
}

#' Helper function to calculate third term of the ELBO: E[log p(z|theta)]
#'
#' @description For given gamma, phi, total cell type counts, total gene counts,
#' and Spatial Transcriptomics data, calculate third term of the ELBO: E[log p(z|theta)]
#'
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param phi.p Variational Inference parameter 3-dimensional phi array,
#'      with dimension pixel x gene x cell type
#' @param ct.count total cell type count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, with dimension pixel x gene
#'
#' @return Calculated third term E[log p(z|theta)]
#'
#' @noRd
#'
#' @keywords internal


t3 <- function(phi.p, gamma.p, pixel.count, ct.count, count.data, gene.count) {
  # E[log p(z|theta)]
  a = digamma(gamma.p) - matrix(rep(digamma(rowSums(gamma.p)), each = ct.count),
                                nrow = pixel.count, byrow = T)
  new_phi = phi.p
  for (i in 1:ct.count) {
    new_phi[,,i] <- as.matrix(new_phi[,,i]*count.data)
  }
  sum(sapply(1:dim(new_phi)[2], function(n) sum(new_phi[, n, ] * a)))
}



#' Helper function to calculate fourth term of the ELBO: E[log p(w|z,beta)]
#'
#' @description For given tau, phi, total cell type, total pixel, total gene count,
#' calculate fourth term of the ELBO: E[log p(w|z,beta)]
#'
#' @param tau.p Variational Inference parameter tau matrix,
#'      with dimension cell type x gene
#' @param phi.p Variational Inference parameter 3-dimensional phi array,
#'      with dimension pixel x gene x cell type
#' @param ct.count total cell type count
#' @param pixel.count total pixel count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, with dimension pixel x gene
#'
#' @return Calculated fourth term E[log p(w|z,beta)]
#'
#' @noRd
#'
#' @keywords internal


t4 <- function(gene.count, pixel.count, ct.count, phi.p, tau.p, count.data) {
  # E[log p(w|z,beta)]
  a = t(digamma(tau.p) - matrix(rep(digamma(rowSums(tau.p)), each = gene.count),
                                nrow = ct.count, byrow = T))
  new_phi = phi.p
  for (i in 1:ct.count) {
    new_phi[,,i] <- as.matrix(new_phi[,,i]*count.data)
  }
  collapsed_matrix <- apply(new_phi, c(2, 3), sum)
  sum(collapsed_matrix * a)
}




#' Helper function to calculate fifth term of the ELBO: E[log q(theta)]
#'
#' @description For given gamma and total pixel count,
#' calculate fifth term of the ELBO: E[log q(theta)]
#'
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param pixel.count total pixel count
#'
#' @return Calculated fifth term E[log q(theta)]
#'
#' @noRd
#'
#' @keywords internal

t5 <- function(pixel.count, gamma.p) {
  # E[log q(theta)]
  val <- 0
  for (i in 1:pixel.count) {
    val <- val + lgamma(sum(gamma.p[i,])) - sum(lgamma(gamma.p[i,]))
    gamma.vec1 <- gamma.p[i,] - 1
    gamma.vec2 <- digamma(gamma.p[i,]) - digamma(sum(gamma.p[i,]))
    val <- val + gamma.vec1 %*% gamma.vec2
  }
  val
}

#' Helper function to calculate sixth term of the ELBO: E[log q(z)]
#'
#' @description For given phi, total pixel count, total cell type count, total
#' gene count and Spatial Transcriptomics data, calculate sixth term of the
#' ELBO: E[log q(z)]
#'
#' @param phi.p Variational Inference parameter 3-dimensional phi array,
#'      with dimension pixel x gene x cell type
#' @param ct.count total cell type count
#' @param pixel.count total pixel count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, with dimension pixel x gene
#'
#' @return Calculated sixth term E[log q(z)]
#'
#' @noRd
#'
#' @keywords internal


t6 <- function(pixel.count, ct.count, gene.count, phi.p, count.data) {
  # E[log q(z)]
  new_phi = phi.p * log(phi.p)
  new_phi[is.nan(new_phi)] = 0
  total_sum = 0
  for (i in 1:ct.count) {
    particular_phi = new_phi[,,i]
    total_sum = total_sum + sum(particular_phi*count.data)

  }
  total_sum
}

#' Helper function to calculate seventh term of the ELBO: E[log q(beta)]
#'
#' @description For tau and total cell type count, calculate seventh term of the
#' ELBO: E[log q(beta)]
#'
#' @param tau.p Variational Inference parameter tau matrix,
#'      with dimension cell type x gene
#' @param ct.count total cell type count
#'
#' @return Calculated seventh term E[log q(beta)]
#'
#' @noRd
#'
#' @keywords internal

t7 <- function(ct.count, tau.p) {
  # E[log q(beta)]
  val <- 0
  for (i in 1:ct.count) {
    val <- lgamma(sum(tau.p[i,])) - sum(lgamma(tau.p[i,]))
    tau.vec1 <- tau.p[i,] - 1
    tau.vec2 <- digamma(tau.p[i,]) - digamma(sum(tau.p[i,]))
    val <- val + tau.vec1 %*% tau.vec2
  }
  val
}



#' Helper function to update parameter phi in the VI algorithm
#'
#' @description For given parameters gamma, tau, total pixel count,
#' total cell type count and total gene count, update the phi parameter
#'
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param tau.p Variational Inference parameter tau matrix,
#'      with dimension cell type x gene
#' @param alpha.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param phi.p Variational Inference parameter 3-dimensional phi array,
#'      with dimension pixel x gene x cell type
#' @param ct.count total cell type count
#' @param pixel.count total pixel count
#' @param gene.count total gene count
#'
#' @return Updated phi parameter, dimension = pixel x gene x cell type
#'
#' @noRd
#'
#' @keywords internal


phi.update <- function(pixel.count, ct.count, gene.count, gamma.p, tau.p) {
  # Better version
  new.phi <- array(0, dim = c(pixel.count, gene.count, ct.count))

  gamma.vec <- digamma(rowSums(gamma.p))
  gamma.mat <- matrix(rep(gamma.vec, each = ct.count), ncol = ct.count,
                      byrow = TRUE)
  new.gamma <- digamma(gamma.p) - gamma.mat

  tau.vec <- digamma(rowSums(tau.p))
  tau.mat <- matrix(rep(tau.vec, each = gene.count), ncol = gene.count,
                    byrow = TRUE)
  new.tau <- digamma(tau.p) - tau.mat

  for (i in 1:pixel.count) {
    gamma.vec <- new.gamma[i,]
    gamma.matrix <- matrix(rep(gamma.vec, each = gene.count), nrow=gene.count)
    tau.matrix <- t(new.tau)
    new.phi[i,,] <- as.matrix(exp(gamma.matrix + tau.matrix))
    new.phi[i,,] <- new.phi[i,,]/rowSums(new.phi[i,,])
  }


  new.phi
}


#' Helper function to update parameter gamma in the VI algorithm
#'
#' @description For given parameters alpha, phi and total pixel count, total
#' gene count, total cell type count and the Spatial Transcriptomics data,
#' update the gamma parameter
#'
#'
#' @param alpha.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param phi.p Variational Inference parameter 3-dimensional phi array,
#'      with dimension pixel x gene x cell type
#' @param ct.count total cell type count
#' @param pixel.count total pixel count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, with dimension pixel x gene
#'
#' @return Updated gamma parameter, dimension = pixel x cell type
#'
#' @noRd
#'
#' @keywords internal

gamma.update <- function(pixel.count, gene.count, ct.count, alpha.p, phi.p,
                         count.data) {
  new.gamma <- matrix(rep(alpha.p, each = pixel.count), nrow = pixel.count,
                      byrow=TRUE)
  for (i in 1:ct.count) {
    new.gamma[,i] <- new.gamma[,i] + rowSums(count.data * phi.p[,,i])
  }
  new.gamma
}


#' Helper function to update parameter tau in the VI algorithm
#'
#' @description For given reference, parameter phi and total cell type count, total
#' gene count, total pixel count and the Spatial Transcriptomics data,
#' update the tau parameter
#'
#' @param ref.p reference matrix, with dimension cell type x gene
#' @param phi.p Variational Inference parameter 3-dimensional phi array,
#'      with dimension pixel x gene x cell type
#' @param ct.count total cell type count
#' @param pixel.count total pixel count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, with dimension pixel x gene
#'
#' @return Updated tau parameter, dimension = cell type x gene
#'
#' @noRd
#'
#' @keywords internal

tau.update <- function(ct.count, gene.count, pixel.count, ref.p, phi.p, count.data) {
  new.tau <- ref.p
  for (i in 1:ct.count) {
    new.tau[i,] <- new.tau[i,] + colSums(count.data * phi.p[,,i])
  }
  new.tau
}

#' Helper function to update parameter alpha in the VI algorithm
#'
#' @description For given parameters alpha and gamma, total pixel count,
#' total gene count, update the alpha parameter
#'
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param alpha.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param ct.count total cell type count
#' @param pixel.count total pixel count
#'
#' @return Updated alpha parameter, length = total cell type count
#'
#' @noRd
#'
#' @keywords internal

alpha.update <- function(pixel.count, ct.count, alpha.p, gamma.p) {
  threshold <- 0.0001
  curr.alpha <- alpha.step(pixel.count, ct.count, alpha.p, gamma.p)
  while(TRUE) {
    new.alpha <- alpha.step(pixel.count, ct.count, curr.alpha, gamma.p)
    if (t2(new.alpha, gamma.p, pixel.count) -
        t2(curr.alpha, gamma.p, pixel.count) < threshold) {
      break
    }
    curr.alpha <- new.alpha
  }
  new.alpha
}



#' Helper function to perform an iterative step during update of alpha
#'
#' @description For given parameters alpha, gamma and total pixel count and total
#' cell type count, perform an iterative update step for alpha
#'
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param alpha.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param ct.count total cell type count
#' @param pixel.count total pixel count
#'
#' @return Updated alpha parameter in an iteration, length = total cell type count
#'
#' @noRd
#'
#' @keywords internal

alpha.step <- function(pixel.count, ct.count, alpha.p, gamma.p) {
  g.vec <- pixel.count * (digamma(sum(alpha.p)) - digamma(alpha.p))

  for (i in 1:pixel.count) {
    for (j in 1:ct.count) {
      g.vec[j] <- g.vec[j] + digamma(gamma.p[i,j]) - digamma(sum(gamma.p[i,]))
    }
  }
  z <- trigamma(sum(alpha.p)) * pixel.count
  h.vec <- pixel.count * trigamma(alpha.p) * (-1)
  c.const <- sum(g.vec/h.vec)/(1/z + sum(1/h.vec))
  update.vec <- (g.vec - c.const)/h.vec
  new.alpha <- alpha.p - update.vec
  const <- 1/2
  while(any(new.alpha <= 0)) {
    new.alpha <- alpha.p - update.vec * const
    const <- const/2
  }
  new.alpha
}







#' Run FlexiDeconv algorithm for cell type deconvolution
#'
#' @description Given the Spatial Transcriptomics data (pixel x gene) and
#' reference (cell type x gene), as well as prior constant for the reference, run
#' FlexiDeconv algorithm to obtain:
#' 1. Cell type proportion per pixel
#' 2. Estimated gene expression profile for each cell type
#'
#' @param spatial pixel x gene matrix representing Spatial Transcriptomics data
#' @param reference cell type x gene matrix representing reference, row sum should be 1
#' @param prior_const vector of length = total number of provided cell type in the reference,
#'      representing weight assigned to each cell type.
#' @param verbose (default: TRUE)
#' @param iter_print_freq (default: 50)
#'
#' @return A list containing estimated cell type proportion per pixel and estimated
#' gene expression profile for each cell type
#'
#' @examples
#'
#' data(mouseHypothalamus)
#' new_reference <- appendPlaceholder(mouseHypothalamus$reference, numPlaceholder=1)
#' prior_const = c(rep(5, 8), 0.1)
#' output = runFlexiDeconv(as.matrix(mouseHypothalamus$spatial),
#' new_reference,
#' prior_const)
#'
#' @export




runFlexiDeconv <- function(spatial, reference, prior_const,
                            verbose = T, iter_print_freq = 50) {

  ## retrieve important statistics about the data
  pixel_num = dim(spatial)[1]
  gene_num = dim(spatial)[2]
  cell_types = rownames(reference)
  medium_molecule_count = median(rowSums(spatial))

  prior_const_matrix = matrix(prior_const, nrow = dim(reference)[1],
                              ncol = dim(reference)[2], byrow = FALSE)
  prior_const_scaled = prior_const_matrix * medium_molecule_count * pixel_num/length(cell_types)



  alpha.vi <- bayess::rdirichlet(n = 1, par = rep(5, length(cell_types)))[1,] * 50
  if (verbose) {
    message("Length of vector alpha: ", length(alpha.vi))
  }


  # gamma is pixel x cell type
  gamma <- length(cell_types) * bayess::rdirichlet(n = pixel_num, par = rep(30, length(cell_types)))
  if (verbose) {
    message("Parameter gamma: ", nrow(gamma), " pixels x ", ncol(gamma), " cell types")
  }

  # phi is pixel x gene x cell type
  phi <- array(0, dim = c(pixel_num, gene_num, length(cell_types)))
  for (i in 1:pixel_num) {
    phi[i,,] = bayess::rdirichlet(n = gene_num, par = rep(30, length(cell_types)))
  }
  if (verbose) {
    message("Parameter phi: ",dim(phi)[1], " pixels x ",dim(phi)[2], " genes x ",
          dim(phi)[3], " cell types")
  }

  # tau is cell type x gene
  tau <- length(cell_types) * bayess::rdirichlet(n = length(cell_types), par = rep(10, gene_num))
  if (verbose) {
    message("Parameter tau: ", nrow(tau), " cell types x ", ncol(tau), " genes")
  }

  ELBOs <- c(ELBO(reference*prior_const_scaled, gamma, tau, alpha.vi, phi, length(cell_types),
                  pixel_num, gene_num, spatial))
  iter <- 1
  while (TRUE) {

    phi <- phi.update(pixel_num, length(cell_types), gene_num, gamma, tau)
    gamma <- gamma.update(pixel_num, gene_num, length(cell_types), alpha.vi, phi,
                          spatial)
    tau <- tau.update(length(cell_types), gene_num, pixel_num,
                      reference*prior_const_scaled, phi, spatial)
    alpha.vi <- alpha.update(pixel_num, length(cell_types), alpha.vi, gamma)
    iter <- iter + 1

    if (iter%%iter_print_freq != 0) {
      next
    }

    ELBOs <- c(ELBOs, ELBO(prior_const_scaled*reference, gamma, tau, alpha.vi, phi, length(cell_types),
                           pixel_num, gene_num, spatial))
    if (verbose) {
      message("Iteration: ", iter,". Curr ELBO: ",ELBOs[length(ELBOs)])
    }
    if (ELBOs[length(ELBOs)] - ELBOs[length(ELBOs) - 1] <= 0.0001) {
      break
    }
  }

  rownames(gamma) = rownames(spatial)
  colnames(gamma) = rownames(reference)
  dimnames(phi) <- list(
    Pixel = rownames(spatial),  # Label for the first dimension (pixels)
    Gene = colnames(spatial),    # Label for the second dimension (genes)
    CellType = rownames(reference)               # Label for the third dimension (cell types)
  )
  rownames(tau) = rownames(reference)
  colnames(tau) = colnames(reference)

  output = list(ELBO = ELBOs, reference = reference, alpha = alpha.vi, gamma = gamma,
                phi = phi, tau = tau, total_iter = iter)
  return(output)
}




















