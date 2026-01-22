ELBO <- function(eta.p, gamma.p, tau.p, alpha.p, phi.p, ct.count, f.p,
                 pixel.count, gene.count, count.data, isoform.dict) {
  start_time <- proc.time()
  term1 <- t1(eta.p, ct.count, tau.p)
  print(paste("1",term1))
  #print("ELBO func: term 1 done")
  term2 <- t2(alpha.p, gamma.p, pixel.count)
  print(paste("2",term2))
  #print("ELBO func: term 2 done")
  term3 <- t3(phi.p, gamma.p, pixel.count, ct.count, count.data, gene.count)
  print(paste("3",term3))
  #print("ELBO func: term 3 done")
  term4 <- t4(gene.count, pixel.count, ct.count, phi.p, tau.p, count.data, isoform.dict)
  print(term4)
  #print("ELBO func: term 4 done")
  term5 <- t5(pixel.count, gamma.p)
  print(term5)
  #print("ELBO func: term 5 done")
  term6 <- t6(pixel.count, ct.count, gene.count, phi.p, count.data)
  print(term6)
  #print("ELBO func: term 6 done")
  term7 <- t7(ct.count, tau.p)
  print(term7)
  #print("ELBO func: term 7 done")
  term8 <- t8(pixel.count, f.p, phi.p, count.data,ct.count)
  print(term8)
  
  gc()
  print(paste("ELBO took:", (proc.time() - start_time)[3], "seconds"))
  term1 + term2 + term3 + term4 - term5 - term6 - term7 + term8
}

t1 <- function(eta.p, ct.count, tau.p) {
  # E[p(beta)]
  val <- 0
  for (i in 1:ct.count) {
    val <- val + lgamma(sum(eta.p[i,])) - sum(lgamma(eta.p[i,]))
    tau.row <- tau.p[i,]
    tau.vec <- digamma(tau.row) - digamma(sum(tau.row))
    val <- val + (eta.p[i,] - 1) %*% tau.vec
  }
  val
}

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

t3 <- function(phi.p, gamma.p, pixel.count, ct.count, count.data, gene.count) {
  # E[log p(z|theta)]
  a = digamma(gamma.p) - matrix(rep(digamma(rowSums(gamma.p)), each = ct.count), 
                                nrow = pixel.count, byrow = T)
  new_phi = phi.p
  for (i in 1:ct.count) {
    #new_phi[,,i] <- as.matrix(new_phi[,,i]*st_data)
    new_phi[,,i] <- as.matrix(new_phi[,,i]*count.data)
  }
  sum(sapply(1:dim(new_phi)[2], function(n) sum(new_phi[, n, ] * a)))
}

t4 <- function(gene.count, pixel.count, ct.count, phi.p, tau.p, count.data, isoform.dict) {
  # E[log p(w|z,beta)]
  
  # gene x cell type
  a = t(digamma(tau.p) - matrix(rep(digamma(rowSums(tau.p)), each = gene.count), 
                                nrow = ct.count, byrow = T))
  new_phi = phi.p
  for (i in 1:ct.count) {
    #new_phi[,,i] <- as.matrix(new_phi[,,i]*st_data)
    new_phi[,,i] <- as.matrix(new_phi[,,i]*count.data)
  }
  # transcript x cell type
  collapsed_matrix <- apply(new_phi, c(2, 3), sum)
  new_collapsed_matrix = matrix(0, nrow = gene.count, ncol = ct.count)
  curr_at = 1
  for (i in 1:gene.count) {
    num_isoform = length(isoform.dict[[i]])
    new_collapsed_matrix[i,] = colSums(collapsed_matrix[curr_at:(curr_at+num_isoform-1),,drop=F])
    curr_at = curr_at + num_isoform
  }
  sum(new_collapsed_matrix * a)
}

#t4(num_gene, num_pixel, num_celltype, phi.vi, tau.vi, spatial, isoform_dict)


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

t6 <- function(pixel.count, ct.count, gene.count, phi.p, count.data) {
  # E[log q(z)]
  new_phi = phi.p * log(phi.p)
  new_phi[is.nan(new_phi)] = 0
  total_sum = 0
  for (i in 1:ct.count) {
    total_sum = total_sum + sum(new_phi[,,i]*count.data)  
  }
  total_sum
}

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

t8 <- function(pixel.count, f.p, phi.p, count.data,ct.count) {
  # E[log p(i|z,w,f)]
  new_phi = phi.p
  f.p.transformed = log(t(f.p))
  val = 0
  for (i in 1:pixel.count) {
    new_phi[i,,] <- new_phi[i,,] * f.p.transformed
  }
  for (k in 1:ct.count) {
    val = val + sum(as.matrix(new_phi[,,k] * count.data))
    print(paste(k,sum(as.matrix(new_phi[,,k] * count.data))))
  }
  val
}

# 
# 
# t4_test <- function(gene.count, pixel.count, ct.count, phi.p, tau.p, count.data, isoform.dict) {
#   val = 0
#   for (i in 1:pixel.count) {
#     for (j in 1:ct.count) {
#       for (k in 1:gene.count) {
#         tau_kn = digamma(tau.p[j,k])
#         tau_sum = digamma(sum(tau.p[j,]))
#         for (l in isoform.dict[[k]]) {
#           spatial_count = count.data[i, paste("g",k,l,sep="")]
#           phi_term = phi.p[i,which(colnames(count.data) == paste("g",k,l,sep="")),j]
#           val = val + spatial_count * phi_term * (tau_kn - tau_sum)
#         }
#       }
#     }
#   }
#   val
# }
# 
# t4_test(num_gene, num_pixel, num_celltype, phi.vi, tau.vi, spatial, isoform_dict)

# this should be the new f.update
log.count = log(spatial)
new.log.phi = log.phi.vi
for (i in 1:34) {
  new.log.phi[,,i] = log.phi.vi[,,i] + as.matrix(log.count)
}
c_val = apply(new.log.phi, MARGIN=c(2,3), FUN=max)
c_3d <- aperm(array(rep(c_val, times = 918), dim = c(1069,34,918)), c(3, 1, 2))
exponential_log = exp(new.log.phi - c_3d)
sum_along_pixel = log(apply(exponential_log, MARGIN=c(2,3), FUN=sum))
new_f_unnormalized = t(c_val + sum_along_pixel)
curr_at = 1
for (i in 1:1000) {
  num_isoform = length(isoform_dict[[i]])
  if (num_isoform == 1) {
    normalizing_factor = matrix(new_f_unnormalized[,curr_at], nrow = 34)
  } else {
    mat_interest = new_f_unnormalized[,curr_at:(curr_at+num_isoform-1)]
    c = apply(mat_interest, 1, max)
    u = exp(mat_interest - matrix(rep(c,num_isoform), nrow = 34, ncol = num_isoform))
    normalizing_factor = matrix(rep(c + log(rowSums(u)), num_isoform), nrow = 34)
  }
  if (curr_at == 1) {
    second_factor = normalizing_factor
  } else {
    second_factor = cbind(second_factor, normalizing_factor)
  }
  curr_at = curr_at + num_isoform
}
f_normalized = new_f_unnormalized - second_factor
#######

f.update <- function(ct.count, isoform.dict, gene.count,
                     phi.p, count.data) {
  # Better version
  new.phi <- phi.p
  for (i in 1:ct.count) {
    new.phi[,,i] = as.matrix(phi.p[,,i] * count.data)
  }
  new.f = t(apply(new.phi, MARGIN = c(2, 3), sum))
  curr_at = 1
  for (i in 1:gene.count) {
    num_isoform = length(isoform.dict[[i]])
    new.f[,curr_at:(curr_at+num_isoform-1)] = new.f[,curr_at:(curr_at+num_isoform-1)]/rowSums(new.f[,curr_at:(curr_at+num_isoform-1),drop=F])
    curr_at = curr_at + num_isoform
  }
  new.f
}




phi.update <- function(pixel.count, ct.count, gene.count, gamma.p, tau.p, isoform.dict, f.p,
                       transcript.count) {
  # Better version
  new.phi <- array(0, dim = c(pixel.count, transcript.count, ct.count))
  
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
    intermediate.new.phi[i,,] <- as.matrix(exp(gamma.matrix + tau.matrix))
    #intermediate.new.phi[i,,] <- intermediate.new.phi[i,,]/rowSums(intermediate.new.phi[i,,])
    #new.phi[i,,] <- new.phi[i,,]/rowSums(new.phi[i,,])
  }
  curr_at = 1
  for (i in 1:gene.count) {
    num_isoform = length(isoform.dict[[i]])
    temp = array(rep(intermediate.new.phi[,i,], each = num_isoform), dim = c(num_isoform, pixel.count, ct.count))
    new.phi[,curr_at:(curr_at+num_isoform-1),] = aperm(temp, c(2, 1, 3))
    #curr_isoform_matrix = t(f.p[,curr_at:(curr_at+num_isoform-1)])
    #temp2 = array(rep(t(curr_isoform_matrix,each=pixel.count)), dim = c(pixel.count, num_isoform, ct.count))
    #new.phi[,curr_at:(curr_at+num_isoform-1),] = new.phi[,curr_at:(curr_at+num_isoform-1),] * curr_isoform_matrix
    curr_at = curr_at + num_isoform
  }
  
  for (i in 1:pixel.count) {
    new.phi[i,,] <- new.phi[i,,] * t(f.p)
    new.phi[i,,] <- new.phi[i,,]/rowSums(new.phi[i,,])
  }
  new.phi
}
#phi.update(256, 6, 120, gamma.vi, tau.vi, isoform_dict, f.vi, 199)


gamma.update <- function(pixel.count, ct.count, alpha.p, phi.p,
                         count.data) {
  new.gamma <- matrix(rep(alpha.p, each = pixel.count), nrow = pixel.count, 
                      byrow=TRUE)
  for (i in 1:ct.count) {
    new.gamma[,i] <- new.gamma[,i] + rowSums(count.data * phi.p[,,i])
  }
  new.gamma
}
#gamma.update(256, 6, alpha.vi, phi.vi, spatial)

tau.update <- function(ct.count, pixel.count, eta.p, phi.p, count.data, transcript.count,
                       isoform.dict, gene.count) {
  # CHECK TAU UPDATE
  intermediate.tau <- matrix(0, nrow = ct.count, ncol = transcript.count)
  new.tau <- eta.p
  for (i in 1:ct.count) {
    intermediate.tau[i,] <- intermediate.tau[i,] + colSums(count.data * phi.p[,,i])
  }
  curr_at = 1
  for (i in 1:gene.count) {
    num_isoform = length(isoform.dict[[i]])
    new.tau[,i] = new.tau[,i] + rowSums(intermediate.tau[,curr_at:(curr_at+num_isoform-1),drop=F])
    curr_at = curr_at + num_isoform
  }
  
  new.tau
}
#tau.update(6, 256, eta, phi.vi, spatial, 199, isoform_dict, 120)


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

alpha.step <- function(pixel.count, ct.count, alpha.p, gamma.p) {
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
  new.alpha
}




logsum <- function(v) {
  # import v is a vector where each entry is of form log(p_i)
  # goal is to compute (log(exp(v1) + ... + exp(vN)))
  c = max(v)
  u = exp(v - c)
  return(log(sum(u)) + c)
}



