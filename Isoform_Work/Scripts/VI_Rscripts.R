
# Now variational inference starts here

ELBO <- function(eta.p, gamma.p, tau.p, alpha.p, phi.p, ct.count,
                 pixel.count, gene.count, count.data) {
  start_time <- proc.time()
  term1 <- t1(eta.p, ct.count, tau.p)
  #print(term1)
  #print("ELBO func: term 1 done")
  term2 <- t2(alpha.p, gamma.p, pixel.count)
  #print(term2)
  #print("ELBO func: term 2 done")
  term3 <- t3(phi.p, gamma.p, pixel.count, ct.count, count.data, gene.count)
  #print(term3)
  #print("ELBO func: term 3 done")
  term4 <- t4(gene.count, pixel.count, ct.count, phi.p, tau.p, count.data)
  #print(term4)
  #print("ELBO func: term 4 done")
  term5 <- t5(pixel.count, gamma.p)
  #print(term5)
  #print("ELBO func: term 5 done")
  term6 <- t6(pixel.count, ct.count, gene.count, phi.p, count.data)
  #print(term6)
  #print("ELBO func: term 6 done")
  term7 <- t7(ct.count, tau.p)
  #print(term7)
  #print("ELBO func: term 7 done")
  gc()
  print(paste("ELBO took:", (proc.time() - start_time)[3], "seconds"))
  term1 + term2 + term3 + term4 - term5 - term6 - term7
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

# t3 <- function(phi.p, gamma.p, pixel.count, ct.count, count.data, gene.count) {
#   val <- 0
#   for (i in 1:pixel.count) {
#     const <- digamma(sum(gamma.p[i,]))
#     for (j in 1:ct.count) {
#       gamma.const <- digamma(gamma.p[i,j]) - const
#
#       local.val <- 0
#       for (k in 1:gene.count) {
#         local.val <- local.val + count.data[i,k] * phi.p[i,k,j]
#       }
#       val <- val + local.val * gamma.const
#     }
#   }
#   val
# }

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

t4 <- function(gene.count, pixel.count, ct.count, phi.p, tau.p, count.data) {
  # E[log p(w|z,beta)]
  a = t(digamma(tau.p) - matrix(rep(digamma(rowSums(tau.p)), each = gene.count),
                                nrow = ct.count, byrow = T))
  new_phi = phi.p
  for (i in 1:ct.count) {
    #new_phi[,,i] <- as.matrix(new_phi[,,i]*st_data)
    new_phi[,,i] <- as.matrix(new_phi[,,i]*count.data)
  }
  collapsed_matrix <- apply(new_phi, c(2, 3), sum)
  sum(collapsed_matrix * a)
}

# t4 <- function(gene.count, pixel.count, ct.count, phi.p, tau.p, count.data) {
#   val <- 0
#   for (i in 1:ct.count) {
#     const <- digamma(sum(tau.p[i,]))
#     for (j in 1:pixel.count) {
#       for (k in 1:gene.count) {
#         # No.molecules of this gene is count.data[j,k]
#         tau.const <- digamma(tau.p[i,k]) - const
#         val <- val + count.data[j,k]*phi.p[j,k,i]*tau.const
#       }
#     }
#   }
#   val
# }

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

# t6 <- function(pixel.count, ct.count, gene.count, phi.p, count.data) {
#   val <- 0
#   for (i in 1:pixel.count) {
#     for (j in 1:gene.count) {
#       number.of.genes <- count.data[i,j]
#       for (k in 1:ct.count) {
#         if (phi.p[i,j,k] <= 0) {
#           val <- val
#         } else {
#           val <- val + number.of.genes * phi.p[i,j,k] * log(phi.p[i,j,k])
#         }
#       }
#     }
#   }
#   val
# }

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

t7 <- function(ct.count, tau.p) {
  # E[log q(beta)]
  val <- 0
  for (i in 1:ct.count) {
    val <- val + lgamma(sum(tau.p[i,])) - sum(lgamma(tau.p[i,]))
    tau.vec1 <- tau.p[i,] - 1
    tau.vec2 <- digamma(tau.p[i,]) - digamma(sum(tau.p[i,]))
    val <- val + tau.vec1 %*% tau.vec2
  }
  val
}


# phi.update <- function(pixel.count, ct.count, gene.count, gamma.p, tau.p) {
#   new.phi <- array(0, dim = c(pixel.count, gene.count, ct.count))
#   for (i in 1:pixel.count) {
#     const <- digamma(sum(gamma.p[i,]))
#     for (j in 1:ct.count) {
#       gamma.const <- digamma(gamma.p[i, j]) - const
#       for (k in 1:gene.count) {
#         tau.const <- digamma(tau.p[j,k]) - digamma(sum(tau.p[j,]))
#         new.phi[i,k,j] <- exp(gamma.const + tau.const)
#       }
#     }
#
#     for (k in 1:gene.count) {
#       new.phi[i,k,] <- new.phi[i,k,]/sum(new.phi[i,k,])
#     }
#   }
#   new.phi
# }

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
    new.phi[i,,] <- as.matrix(gamma.matrix + tau.matrix)
    #new.phi[i,,] <- new.phi[i,,]/rowSums(new.phi[i,,])
    #new.phi[i,,] <- new.phi[i,,]
  }
  x <- new.phi
  c <- apply(x, c(1, 2), max)
  extended_c <- array(rep(c, times = length(cell_types)), dim = c(pixel_num, gene_num, length(cell_types)))
  exp_val <- exp(x - extended_c)
  log_sum_up <- log(apply(exp_val, c(1,2), sum))
  extended_log_sum_up <- array(rep(log_sum_up, times = length(cell_types)), dim = c(pixel_num, gene_num, length(cell_types)))
  y = extended_c + extended_log_sum_up
  new_phi = x - y

  new_phi
}

phi.update.old <- function(pixel.count, ct.count, gene.count, gamma.p, tau.p) {
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

# gamma.update <- function(pixel.count, gene.count, ct.count, alpha.p, phi.p,
#                          count.data){
#   new.gamma <- matrix(1, nrow = pixel.count, ncol = ct.count)
#   for (i in 1:pixel.count) {
#     for (j in 1:ct.count) {
#       new.gamma[i,j] <- alpha.p[j] + phi.p[i,,j] %*% count.data[i,]
#     }
#   }
#   new.gamma
# }

gamma.update <- function(pixel.count, gene.count, ct.count, alpha.p, phi.p,
                         count.data) {
  new.gamma <- matrix(rep(alpha.p, each = pixel.count), nrow = pixel.count,
                      byrow=TRUE)
  for (i in 1:ct.count) {
    new.gamma[,i] <- new.gamma[,i] + rowSums(count.data * phi.p[,,i])
  }
  new.gamma
}

# tau.update <- function(ct.count, gene.count, pixel.count, eta.p, phi.p, count.data) {
#   new.tau <- matrix(1, nrow = ct.count, ncol = gene.count)
#   for (i in 1:ct.count) {
#     for (j in 1:gene.count) {
#       const <- 0
#       for (k in 1:pixel.count) {
#         const <- const + count.data[k,j] * phi.p[k,j,i]
#       }
#       new.tau[i,j] <- const + eta.p[i,j]
#     }
#   }
#   new.tau
# }


# tau.update <- function(ct.count, gene.count, pixel.count, eta.p, phi.p, count.data) {
#   # Better version
#   term2.step1 <- phi.p
#   for (k in 1:ct.count) {
#     term2.step1[,,k] <- term2.step1[,,k] * count.data
#   }
#   term2.step2 <- apply(term2.step1, c(2,3), sum)
#   t(term2.step2) + eta.p
# }

tau.update <- function(ct.count, gene.count, pixel.count, eta.p, phi.p, count.data) {
  new.tau <- eta.p
  for (i in 1:ct.count) {
    new.tau[i,] <- new.tau[i,] + colSums(count.data * phi.p[,,i])
  }
  new.tau
}



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


