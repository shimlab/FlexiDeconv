############################ ELBO FUNCTIONS ###################################
#' Helper function to calculate the ELBO
#' 
#' @description For given reference, gamma, tau, alpha, phi parameter values, as well
#' as cell type counts, pixel counts, gene counts, the spatial transcriptomics
#' data and the isoform distribution, calculate the ELBO
#' 
#' @param ref.p reference matrix, with dimension cell type x gene
#' @param alpha.p Dirichlet hyperparameter vector for cell type proportions,
#'      with vector length cell type
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param log.phi.p Variational Inference parameter 3-dimensional log.phi array,
#'      with dimension pixel x isoform x cell type, where isoforms are indexed
#'      globally across all genes (flattened representation)
#' @param tau.p Variational Inference parameter tau matrix,
#'      with dimension cell type x gene
#' @param log.f.p Variational Inference parameter 2-dimensional log.f array,
#'      with dimension cell type x isoform
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#' @param isoform.dict List mapping genes to their isoforms,
#'      with length gene
#' 
#' @return ELBO value for current parameter set
#' 
#' @noRd
#' 
#' @keywords internal

ELBO <- function(ref.p, alpha.p, gamma.p, log.phi.p, tau.p, log.f.p, pixel.count, 
                 ct.count, gene.count, count.data, isoform.dict, verbose = T) {
  start_time <- proc.time()
  term1 <- t1(ref.p, tau.p, ct.count)
  term2 <- t2(alpha.p, gamma.p, pixel.count)
  term3 <- t3(gamma.p, log.phi.p, pixel.count, ct.count, gene.count, count.data)
  term4 <- t4(log.phi.p, tau.p, pixel.count, ct.count, gene.count, count.data, isoform.dict)
  term5 <- t5(gamma.p, pixel.count)
  term6 <- t6(log.phi.p, pixel.count, ct.count, gene.count, count.data)
  term7 <- t7(tau.p, ct.count)
  # E[log p(i | z,w,f)]
  term8 <- t8(log.phi.p, log.f.p, pixel.count, ct.count, count.data)
  
  gc()
  print(paste("ELBO took: ", (proc.time() - start_time)[3], "seconds"))
  
  if (verbose) {
    print(paste("t1: E[p(beta)] | ", term1))
    print(paste("t2: E[p(theta | alpha)] | ", term2))
    print(paste("t3: E[p(z | theta)] | ", term3))
    print(paste("t4: E[p(w | z, beta)] | ", term4))
    print(paste("t8: E[p(i | z, w, f)] | ", term8))
    print(paste("t5: E[q(theta)] | ", term5))
    print(paste("t6: E[q(z)] | ", term6))
    print(paste("t7: E[q(beta)] | ", term7))
  }
  
  return(term1 + term2 + term3 + term4 + term8 - term5 - term6 - term7)
}

#### Term 1 ####
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

t1 <- function(ref.p, tau.p, ct.count) {
  # E[p(beta)]
  val <- 0
  
  # Sums over all cell types (i.e. sum_{k = 1}^K)
  for (i in 1:ct.count) {
    # Computes: log Gamma(sum_{n=1}^N b_k eta_{kn}) - sum_{n=1}^N log Gamma(b_k eta_{kn})
    val <- val + lgamma(sum(ref.p[i,])) - sum(lgamma(ref.p[i,]))
    tau.row <- tau.p[i,]
    tau.vec <- digamma(tau.row) - digamma(sum(tau.row))
    # Computes: sum_{n=1}^N [(b_k eta_{kn} - 1) * (digamma(tau_{kn}) - digamma(sum_{j=1}^N tau_{kj}))]
    val <- val + (ref.p[i,] - 1) %*% tau.vec
  }
  return(val)
}


#### Term 2 ####
#' Helper function to calculate second term of the ELBO: E[log p(theta|alpha)]
#'
#' @description For given gamma, alpha, total pixel counts, calculate second
#' term of the ELBO: E[log p(theta|alpha)]
#'
#' @param alpha.p Dirichlet hyperparameter vector for cell type proportions,
#'      with vector length cell type
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param pixel.count total pixel count
#'
#' @return Calculated second term E[log p(theta|alpha)]
#'
#' @noRd
#'
#' @keywords internal

t2 <- function(alpha.p, gamma.p, pixel.count){
  # E[log p(theta|alpha)]
  
  # Computes: sum_{d = 1}^D [log Gamma(sum_{k=1}^K alpha_k) - sum_{k=1}^K log Gamma(alpha_k})]
  val <- pixel.count * (lgamma(sum(alpha.p)) - sum(lgamma(alpha.p)))
  
  
  # Sums over all pixels (i.e. sum_{d = 1}^D)
  for (i in 1:pixel.count) {
    gamma.row <- gamma.p[i,]
    gamma.vec <- digamma(gamma.row) - digamma(sum(gamma.row))
    # Computes: sum_{k=1}^K [(alpha_k - 1) * (digamma(gamma_{dk}) - digamma(sum_{j=1}^K tau_{dj}))]
    val <- val + (alpha.p - 1) %*% gamma.vec
  }
  return(val)
}


#### Term 3 ####
#' Helper function to calculate third term of the ELBO: E[log p(z|theta)]
#'
#' @description For given gamma, phi, total cell type counts, total gene counts,
#' and Spatial Transcriptomics data, calculate third term of the ELBO: E[log p(z|theta)]
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param log.phi.p Variational Inference parameter 3-dimensional log.phi array,
#'      with dimension pixel x isoform x cell type, where isoforms are indexed
#'      globally across all genes (flattened representation)
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#'
#' @return Calculated third term E[log p(z|theta)]
#'
#' @noRd
#'
#' @keywords internal

t3 <- function(gamma.p, log.phi.p, pixel.count, ct.count, gene.count, count.data){
  # Computes digamma(gamma_{dk}) - digamma(sum_{j = 1}^K gamma_{dj})
  # Output: Matrix with dimension pixel.count (D) x ct.count (K)
  a = digamma(gamma.p) - matrix(rep(digamma(rowSums(gamma.p)), each = ct.count),
                                nrow = pixel.count, byrow = T)
  
  new.phi = exp(log.phi.p)
  # Sums over all cell types (i.e. sum_{k = 1}^K)
  for (i in 1:ct.count){
    # Computes: phi_{d n_s k} * count_{d n_s} for cell type k
    # Output: Matrix with dimension pixel.count (D) x total.isoforms
    new.phi[,,i] <- as.matrix(new.phi[,,i] * count.data)
  }
  
  # Sums across all genes (N) for phi * a * count.data
  # Output: Value
  return(sum(sapply(1:dim(new.phi)[2], function(n) sum(new.phi[, n, ] * a))))
}


#### Term 4 ####
#' Helper function to calculate fourth term of the ELBO: E[log p(w|z,beta)]
#'
#' @description For given tau, phi, total cell type, total pixel, total gene count
#' and isoform dict calculate fourth term of the ELBO: E[log p(w|z,beta)]
#'
#' @param log.phi.p Variational Inference parameter 3-dimensional log.phi array,
#'      with dimension pixel x isoform x cell type, where isoforms are indexed
#'      globally across all genes (flattened representation)
#' @param tau.p Variational Inference parameter tau matrix,
#'      with dimension cell type x gene
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#' @param isoform.dict List mapping genes to their isoforms,
#'      with length gene
#'
#' @return Calculated fourth term E[log p(w|z,beta)]
#'
#' @noRd
#'
#' @keywords internal

t4 <- function(log.phi.p, tau.p, pixel.count, ct.count, gene.count, 
               count.data, isoform.dict){
  # E[log p(w|z,beta)]
  
  # Computes digamma(gamma_{dk}) - digamma(sum_{j = 1}^K gamma_{dj})
  # Output: Matrix with dimension total.isoform x ct.count (K)
  a = t(digamma(tau.p) - matrix(rep(digamma(rowSums(tau.p)), each = gene.count), 
                                nrow = ct.count, byrow = T))
  
  new.phi = exp(log.phi.p)
  # Sums over all cell types (i.e. sum_{k = 1}^K)
  for (i in 1:ct.count) {
    new.phi[,,i] <- as.matrix(new.phi[,,i] * count.data)
  }
  
  ########### CHECK isoform.dict structure ###########
  isoform.count.expected = sum(sapply(isoform.dict, length))
  isoform.count.actual = dim(log.phi.p)[2]
  if (isoform.count.expected != isoform.count.actual){
    print(paste("Dimension mismatch! Expected: ", isoform.count.expected,
                "| Actual: ", isoform.count.actual))
  }
  ########### END CHECK ###########
  
  # Collapses pixel dimension in new.phi
  # Output: Matrix with dimension gene.count (N) x ct.count (K) 
  collapsed_matrix <- apply(new.phi, c(2, 3), sum)
  new_collapsed_matrix = matrix(0, nrow = gene.count, ncol = ct.count)
  curr_at = 1
  
  # Sums over all genes (i.e. sum_{n = 1}^N)
  for (i in 1:gene.count) {
    # Computes phi_{dn_s k} * count_{dn_s}
    # Output: Matrix with dimension pixel.count (D) x transcript (T_n)
    num_isoform = length(isoform.dict[[i]])
    new_collapsed_matrix[i,] = colSums(collapsed_matrix[curr_at:(curr_at+num_isoform-1),,drop=F])
    curr_at = curr_at + num_isoform
  }
  return(sum(new_collapsed_matrix * a))
}


#### Term 5 ####
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

t5 <- function(gamma.p, pixel.count) {
  # E[log q(theta)]
  val <- 0
  
  # Sums over all pixels (i.e. sum_{d = 1}^D)
  for (i in 1:pixel.count) {
    # Computes: log Gamma(sum_{k=1}^K gamma_{dk}) - sum_{k=1}^K log Gamma(gamma_{dk})
    val <- val + lgamma(sum(gamma.p[i,])) - sum(lgamma(gamma.p[i,]))
    gamma.vec1 <- gamma.p[i,] - 1
    gamma.vec2 <- digamma(gamma.p[i,]) - digamma(sum(gamma.p[i,]))
    
    # Computes: sum_{k=1}^K [(gamma_{dk} - 1) * (digamma(gamma_{dk}) - digamma(sum_{j=1}^K gamma_{dj}))]
    val <- val + gamma.vec1 %*% gamma.vec2
  }
  return(val)
}


#### Term 6 ####
#' Helper function to calculate sixth term of the ELBO: E[log q(z)]
#'
#' @description For given phi, total pixel count, total cell type count, total
#' gene count and Spatial Transcriptomics data, calculate sixth term of the
#' ELBO: E[log q(z)]
#'
#' @param log.phi.p Variational Inference parameter 3-dimensional log.phi array,
#'      with dimension pixel x isoform x cell type, where isoforms are indexed
#'      globally across all genes (flattened representation)
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#'
#' @return Calculated sixth term E[log q(z)]
#'
#' @noRd
#'
#' @keywords internal

t6 <- function(log.phi.p, pixel.count, ct.count, gene.count, count.data) {
  # E[log q(z)]
  
  # Computes phi_{d n_s k} * log(phi_{d n_s k})
  new.phi = exp(log.phi.p) * log.phi.p
  new.phi[is.nan(new.phi)] = 0
  total_sum = 0
  
  # Sums over all cell types (i.e. sum_{k = 1}^K)
  for (i in 1:ct.count) {
    # Computes: phi_{d n_s k} * log(phi_{d n_s k} * count_{d n_s}
    total_sum = total_sum + sum(new.phi[,,i] * count.data)  
  }
  return(total_sum)
}


#### Term 7 ####
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

t7 <- function(tau.p, ct.count) {
  # E[log q(beta)]
  val <- 0
  
  # Sums over all cell types (i.e. sum_{k = 1}^K)
  for (i in 1:ct.count) {
    # Computes: log Gamma(sum_{n=1}^N tau_{kn}) - sum_{n=1}^N log Gamma(tau_{kn})
    val <- val + lgamma(sum(tau.p[i,])) - sum(lgamma(tau.p[i,]))
    tau.vec1 <- tau.p[i,] - 1
    tau.vec2 <- digamma(tau.p[i,]) - digamma(sum(tau.p[i,]))
    # Computes: sum_{n=1}^N [(tau_{kn} - 1) * (digamma(tau_{kn}) - digamma(sum_{j=1}^N tau_{kj}))]
    val <- val + tau.vec1 %*% tau.vec2
  }
  return(val)
}

#### Term 8 ####
#' Helper function to calculate eighth term of the ELBO: E[log p(i|z,w,f)]
#'
#' @description For given phi, f, total pixel count, total cell type count, 
#'      and Spatial Transcriptomics data, calculate eighth term of the
#'      ELBO: E[log p(i|z,w,f)]
#'
#' @param log.phi.p Variational Inference parameter 3-dimensional log.phi array,
#'      with dimension pixel x isoform x cell type, where isoforms are indexed
#'      globally across all genes (flattened representation)
#' @param log.f.p Variational Inference parameter 2-dimensional log.f array,
#'      with dimension cell type x isoform
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#'
#' @return Calculated sixth term E[log q(z)]
#'
#' @noRd
#'
#' @keywords internal


t8 <- function(log.phi.p, log.f.p, pixel.count, ct.count, count.data) {
  # E[log p(i|z,w,f)]
  new.phi = exp(log.phi.p)
  
  f.p.transformed = t(log.f.p)
  val = 0
  
  # ============================================
  # VALIDATION: Check log.f.p is valid
  # ============================================
  if (any(is.nan(log.f.p))) {
    cat("✗ ERROR in t8: log.f.p contains NaN values!\n")
    cat("  NaN count:", sum(is.nan(log.f.p)), "out of", length(log.f.p), "\n")
    return(NaN)
  }
  
  if (any(log.f.p == Inf)) {
    cat("✗ ERROR in t8: log.f.p contains Inf values!\n")
    cat("  Inf count:", sum(is.infinite(log.f.p)), "\n")
    return(NaN)
  }
   
  # ============================================
  
  # Sums over all pixels (i.e. sum_{d = 1}^D)
  for (i in 1:pixel.count) {
    # Computes: phi_{d n_s k} * log_kn (n_s)
    new.phi[i,,] <- new.phi[i,,] * f.p.transformed
  }
  
  # Sums over all cell types (i.e. sum_{k = 1}^K)
  for (k in 1:ct.count) {
    # Computes: phi_{d n_s k} * log_kn (n_s) * count_{d n_s}
    val = val + sum(as.matrix(new.phi[,,k] * count.data))
    # print(paste(k,sum(as.matrix(new.phi[,,k] * count.data))))
  }
  return(val)
}

######################### END ELBO FUNCTIONS ##################################

##########################  HELPER FUNCTIONS  #################################

#### log-sum-exp trick ####
#' Compute log(sum(exp(x))) in numerically stable way
#'
#' @description Uses the log-sum-exp trick to prevent overflow/underflow when
#'      computing log of sum of exponentials. Filters out non-finite values.
#'
#' @param x Numeric vector of log-space probabilities, 
#'      where x_i = log(p_i)
#'
#' @return Scalar value: log(sum(exp(x)))
#'
#' @references https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
#'
#' @noRd
#'
#' @keywords internal

logsumexp <- function(x) {
  # Handling edge cases, non-finite values (Inf, -Inf, NaN) are filtered out
  x = x[is.finite(x)]
  if (length(x) == 0) return(-Inf)
  if (length(x) == 1) return(x)
  
  x_max = max(x)
  output = x_max + log(sum(exp(x - x_max)))
  return(output)
}


#' Compute log(sum(exp(x))) row-wise for a matrix
#'
#' @description Applies log-sum-exp trick to each row of a matrix independently,
#'      where for each row i it computes log(sum(exp(X[i,]))).
#'
#' @param X Numeric matrix with dimensions (n x m), 
#'      where each entry is in log-space
#'
#' @return Numeric vector of length n, where result[i] = logsumexp(X[i, ])
#'
#' @noRd
#'
#' @keywords internal

#### log-sum-exp row-wise ####
logsumexp.rows <- function(X){
  if (!is.matrix(X)){
    stop("Input should be a matrix")
  }
  
  X[!is.finite(X)] = -Inf
  x_max = apply(X, 1, max)
  x_shifted = sweep(X, 1, x_max, "-")
  output = x_max + log(rowSums(exp(x_shifted)))
  
  output[x_max == -Inf] = -Inf
  
  return(output)
}

#### log-sum-exp col-wise ####
#' Compute log(sum(exp(x))) column-wise for a matrix
#'
#' @description Applies log-sum-exp trick to each column of a matrix independently.
#' For each column j, computes log(sum(exp(X[,j]))).
#'
#' @param X Numeric matrix with dimensions (n x m), where each entry is in log-space
#'
#' @return Numeric vector of length m, where result[j] = logsumexp(X[, j])
#'
#' @noRd
#'
#' @keywords internal

logsumexp.cols <- function(X){
  if (!is.matrix(X)){
    stop("Input should be a matrix")
  }
  
  X[!is.finite(X)] = -Inf
  x_max = apply(X, 2, max)
  x_shifted = sweep(X, 2, x_max, "-")
  output = x_max + log(colSums(exp(x_shifted)))
  
  output[x_max == -Inf] = -Inf
  return(output)
}


#### weighted-log-sum-exp ####
#' Compute log(Σ_i weights[i] * exp(x[i])) using log(sum(exp(x)))
#'
#' @description Weighted log-sum-exp for computing log of weighted sums
#' of exponentials. Handles zero weights and -Inf values correctly.
#'
#' @param x Numeric vector of log-values
#' @param weights Numeric vector of weights (must be non-negative)
#'
#' @return Scalar: log(Σ_i weights[i] * exp(x[i]))
#'
#' @noRd
#'
#' @keywords internal

weighted.logsumexp <- function(x, weights) {
  x[!is.finite(x)] = -Inf
  weights[!is.finite(x)] = 0
  
  max.x = max(x)
  weighted.sum <- sum(weights * exp(x - max.x))
  return(max.x + log(weighted.sum))
}

########################  END HELPER FUNCTIONS  ###############################

##########################  PARAMETER UPDATES  ################################

#### f parameter update ####
#' Helper function to update parameter log.f in the VI algorithm
#'
#' @description For given parameters phi, total cell type count, total gene count,
#'      isoform dict distribution and Spatial Transcriptomics data, 
#'      update the f parameter
#'
#' @param log.phi.p Variational Inference parameter 3-dimensional log.phi array,
#'      with dimension pixel x isoform x cell type, where isoforms are indexed
#'      globally across all genes (flattened representation)
#' @param ct.count total cell type count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#' @param isoform.dict List mapping genes to their isoforms,
#'      with length gene
#'      
#' @return Updated log.f parameter matrix, with dimension cell type x isoform
#'
#' @noRd
#'
#' @keywords internal

f.update <- function(log.phi.p, ct.count, gene.count, count.data, isoform.dict) {
  
  ########### CHECK isoform.dict structure ###########
  isoform.count.expected = sum(sapply(isoform.dict, length))
  isoform.count.actual = dim(log.phi.p)[2]
  if (isoform.count.expected != isoform.count.actual){
    print(paste("Dimension mismatch! Expected: ", isoform.count.expected,
                "| Actual: ", isoform.count.actual))
  }
  ########### END CHECK ###########
  
  isoform.count = dim(log.phi.p)[2]
  log.new.f <- matrix(-Inf, nrow = ct.count, ncol = isoform.count)
  
  for (i in 1:ct.count) {
    for (j in 1:isoform.count){
      log.phi.vec = log.phi.p[, j, i]
      counts.vec <- count.data[, j]
      
      # Collapses pixel dimension in new.f
      # Computes log (phi_{d n_s k} * count_{d n_s}) for cell type k
      # Output: Matrix with dimension  ct.count (K) x isoform
      log.new.f[i, j] = weighted.logsumexp(log.phi.vec, counts.vec)
    }
  }
  
  
  curr_at <- 1
  for (n in 1:gene.count) {
    num_isoform <- length(isoform.dict[[n]])
    isoform.cols <- curr_at:(curr_at + num_isoform - 1)
    
    log.gene.counts <- log.new.f[, isoform.cols, drop = FALSE]
    
    # Normalize each row using logsumexp
    for (k in 1:ct.count){
      log.vals = log.gene.counts[k, ]
      
      # Edge case: If all counts are 0
      # Use uniform distribution (no evidence to prefer any isoform)
      if (all(!is.finite(log.vals))){
        log.new.f[k, isoform.cols] = log(1 / num_isoform)
      } else {
        # Normalize in log-space to avoid underflow
        log.normalizer = logsumexp(log.vals)
        log.new.f[k, isoform.cols] = log.vals - log.normalizer 
      }
    }
    
    curr_at = curr_at + num_isoform
  }
  return(log.new.f)
}


#### phi parameter update ####
#' Helper function to update parameter log.phi in the VI algorithm
#'
#' @description For given parameters gamma, tau, total pixel count,
#' total cell type count and total gene count, update the phi parameter
#'
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param tau.p Variational Inference parameter tau matrix,
#'      with dimension cell type x gene
#' @param log.f.p Variational Inference parameter 2-dimensional log.f array,
#'      with dimension cell type x isoform
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#' @param isoform.dict List mapping genes to their isoforms,
#'      with length gene
#'
#' @return Updated log.phi parameter, dimension = pixel x gene x cell type
#'
#' @noRd
#'
#' @keywords internal

phi.update <- function(gamma.p, tau.p, log.f.p, pixel.count, ct.count, gene.count, 
                       count.data, isoform.dict) {
  # Better version
  transcript.count = dim(count.data)[2]
  output.log.phi <- array(0, dim = c(pixel.count, transcript.count, ct.count))
  log.new.phi <- array(0, dim = c(pixel.count, transcript.count, ct.count))
  
  intermediate.new.phi = array(0, dim = c(pixel.count, gene.count, ct.count))
  
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
    
    # Attempted implementation of logsumexp trick
    intermediate.new.phi[i,,] = as.matrix(gamma.matrix + tau.matrix)
  }
  
  curr_at = 1
  for (i in 1:gene.count) {
    num_isoform = length(isoform.dict[[i]])
    temp = array(rep(intermediate.new.phi[,i,], each = num_isoform), dim = c(num_isoform, pixel.count, ct.count))
    log.new.phi[,curr_at:(curr_at+num_isoform-1),] = aperm(temp, c(2, 1, 3))
    curr_at = curr_at + num_isoform
  }
  
  log.f.p = t(log.f.p)
  log.f.p[!is.finite(log.f.p)] = -Inf # Handling NaN values and setting to p = 0 <=> log p = -Inf
  
  for (i in 1:pixel.count) {
    log.phi.weighted = log.new.phi[i,,] + log.f.p # Calculating new.phi[i,, ]
    log.normalizers = logsumexp.rows(log.phi.weighted) # Equivalent to rowSums(new.phi[...]) in original implementation
    
    # Normalise each row (isoform) across cell types using log-sum-exp
    log.phi.normalized = log.phi.weighted - log.normalizers # Log space form of new.phi[i,,]/rowSums(new.phi[i,,])
    output.log.phi[i,,] = log.phi.normalized
  }
  return(output.log.phi)
}

#### gamma parameter update ####
#' Helper function to update parameter gamma in the VI algorithm
#'
#' @description For given parameters alpha, phi and total pixel count, total
#' gene count, total cell type count and the Spatial Transcriptomics data,
#' update the gamma parameter
#'
#'
#' @param alpha.p Dirichlet hyperparameter vector for cell type proportions,
#'      with vector length cell type
#' @param log.phi.p Variational Inference parameter 3-dimensional log.phi array,
#'      with dimension pixel x isoform x cell type, where isoforms are indexed
#'      globally across all genes (flattened representation)
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#'
#' @return Updated gamma parameter, dimension = pixel x cell type
#'
#' @noRd
#'
#' @keywords internal

gamma.update <- function(alpha.p, log.phi.p, pixel.count, ct.count, 
                         count.data) {
  new.gamma <- matrix(rep(alpha.p, each = pixel.count), nrow = pixel.count, 
                      byrow=TRUE)
  for (i in 1:ct.count) {
    new.gamma[,i] <- new.gamma[,i] + rowSums(count.data * exp(log.phi.p[,,i]))
  }
  return(new.gamma)
}


#### tau parameter update ####
#' Helper function to update parameter tau in the VI algorithm
#'
#' @description For given reference, parameter phi and total cell type count, total
#' gene count, total pixel count and the Spatial Transcriptomics data,
#' update the tau parameter
#'
#' @param ref.p reference matrix, with dimension cell type x gene
#' @param log.phi.p Variational Inference parameter 3-dimensional log.phi array,
#'      with dimension pixel x isoform x cell type, where isoforms are indexed
#'      globally across all genes (flattened representation)
#' @param ct.count total cell type count
#' @param gene.count total gene count
#' @param count.data Spatial Transcriptomics data matrix, 
#'      with dimension pixel x isoform
#' @param isoform.dict List mapping genes to their isoforms,
#'      with length gene
#'
#' @return Updated tau parameter, dimension = cell type x gene
#'
#' @noRd
#'
#' @keywords internal

tau.update <- function(ref.p, log.phi.p, ct.count, gene.count, 
                       count.data, isoform.dict) {
  # CHECK TAU UPDATE
  transcript.count = dim(count.data)[2]
  intermediate.tau <- matrix(0, nrow = ct.count, ncol = transcript.count)
  new.tau <- ref.p
  for (i in 1:ct.count) {
    intermediate.tau[i,] <- colSums(count.data * exp(log.phi.p[,,i]))
  }
  curr_at = 1
  for (i in 1:gene.count) {
    num_isoform = length(isoform.dict[[i]])
    new.tau[,i] = new.tau[,i] + rowSums(intermediate.tau[,curr_at:(curr_at+num_isoform-1),drop=F])
    curr_at = curr_at + num_isoform
  }
  return(new.tau)
}


#### alpha parameter update ####
#' Helper function to update parameter alpha in the VI algorithm
#'
#' @description For given parameters alpha and gamma, total pixel count,
#' total gene count, update the alpha parameter
#'
#' @param alpha.p Dirichlet hyperparameter vector for cell type proportions,
#'      with vector length cell type
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#' @param threshold Convergence threshold for Newton-Raphson iteration,
#'      stops when ELBO improvement less than value (Default: 0.0001)
#'
#' @return Updated alpha parameter, length = total cell type count
#'
#' @noRd
#'
#' @keywords internal

alpha.update <- function(alpha.p, gamma.p, pixel.count, ct.count, 
                         threshold = 0.0001) {
  curr.alpha <- alpha.step(alpha.p, gamma.p, pixel.count, ct.count)
  while(TRUE) {
    new.alpha <- alpha.step(curr.alpha, gamma.p, pixel.count, ct.count)
    # Compares ELBO term t2 = E[log p (theta|alpha)] across both values
    if (t2(new.alpha, gamma.p, pixel.count) - 
        t2(curr.alpha, gamma.p, pixel.count) < threshold) {
      break
    }
    curr.alpha <- new.alpha
  }
  return(new.alpha)
}

#### alpha step function ####
#' Helper function to perform an iterative step during update of alpha
#'
#' @param alpha.p Dirichlet hyperparameter vector for cell type proportions,
#'      with vector length cell type
#' @param gamma.p Variational Inference parameter gamma matrix,
#'      with dimension pixel x cell type
#' @param pixel.count total pixel count
#' @param ct.count total cell type count
#'
#' @return Updated alpha parameter in an iteration, length = total cell type count
#'
#' @noRd
#'
#' @keywords internal

alpha.step <- function(alpha.p, gamma.p, pixel.count, ct.count) {
  g.vec <- pixel.count * (digamma(sum(alpha.p)) - digamma(alpha.p))
  
  for (i in 1:pixel.count) {
    val = digamma(sum(gamma.p[i,]))
    for (j in 1:ct.count) {
      g.vec[j] <- g.vec[j] + digamma(gamma.p[i,j]) - val
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
  return(new.alpha)
}


########################  END PARAMETER UPDATES  ##############################
