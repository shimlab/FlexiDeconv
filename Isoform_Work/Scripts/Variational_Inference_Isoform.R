library(bayess)
library(Matrix)
source("/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Scripts/functions(new).R")

args <- commandArgs(trailingOnly = TRUE)

# Lebrigand Isoform:
spatial <- readRDS(args[1])
reference <- readRDS(args[2])
isoform_dict <- readRDS(args[3])
prior_constant_data = readRDS(args[4])
prior_multiplier = as.numeric(args[5])

output_directory = args[6]

# Define parameters
num_pixel = dim(spatial)[1]
num_gene = dim(reference)[2]
num_celltype = dim(reference)[1]
num_transcript = dim(spatial)[2]
cell_types = rownames(reference)
genes = colnames(reference)

isoform_vec <- c()
for (i in 1:num_gene) {
  isoform_vec = c(isoform_vec, paste(genes[i],isoform_dict[[i]],sep=""))
}


# Read Prior Constant and preprocess
medium_molecule_count = median(rowSums(spatial))
prior_constant_data = prior_multiplier * prior_constant_data * medium_molecule_count * num_pixel/num_celltype


################################################################################
source("/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Scripts/VI_Rscripts_Isoform_Logspace.R")
set.seed(as.numeric(args[7]))

# Modification: Import Gene LR to see if starting on correct estimation helps with VI output
# Gene LR
gene_lr = readRDS("/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Outputs/report_results/pw_2.0/gene_LR_pw_2.0_t1.rds")
gene_estimate_gamma = gene_lr$gamma
rownames(gene_estimate_gamma) = NULL
colnames(gene_estimate_gamma) = NULL

alpha.vi <- rdirichlet(n = 1, par = rep(5, num_celltype))[1,] * 50
# gamma.vi <- num_celltype * rdirichlet(n = num_pixel, par = rep(30, num_celltype))
gamma.vi = gene_estimate_gamma

phi.vi <- array(0, dim = c(num_pixel, num_transcript, num_celltype))
for (i in 1:num_pixel) {
  phi.vi[i,,] = rdirichlet(n = num_transcript, par = rep(30, num_celltype))
}

tau.vi <- num_celltype * rdirichlet(n = num_celltype, par = rep(10, num_gene))

#f.vi <-  matrix(vector("list", num_celltype * num_gene), nrow = num_celltype, ncol = num_gene)
#rownames(f.vi) = cell_types
#colnames(f.vi) = genes
#for (i in 1:num_celltype) {
#  for (j in 1:num_gene) {
#    vec = rdirichlet(n = 1, par = rep(1, length(isoform_dict[[paste("g",j,sep="")]])))[1,]
#    names(vec) = isoform_dict[[paste("g",j,sep="")]]
#    f.vi[[i,j]] = vec
#  }
#}
f.vi <- matrix(0, nrow = num_celltype, ncol = num_transcript)
rownames(f.vi) = cell_types
colnames(f.vi) = isoform_vec
for (i in 1:num_celltype) {
  for (j in 1:num_gene) {
    curr_gene_isoform = isoform_dict[[paste(genes[j],sep="")]]
    #print(f.vi[cell_types[i], paste(genes[j],curr_gene_isoform, sep="")])
    #print(isoform_dict[[paste("g",j,sep="")]])
    f.vi[cell_types[i], paste(genes[j],curr_gene_isoform, sep="")] = rdirichlet(n = 1, 
                                                                                par = rep(1, length(isoform_dict[[genes[j]]])))[1,]
  }
}
log.f.vi = log(f.vi)
log.f.vi[!is.finite(log.f.vi)] = -Inf
log.phi.vi = log(phi.vi)
log.phi.vi[!is.finite(log.phi.vi)] = -Inf
eta <- reference


ELBOs <- c(ELBO(eta*prior_constant_data, alpha.vi, gamma.vi, log.phi.vi,  tau.vi, 
                log.f.vi, num_pixel, num_celltype, num_gene, spatial, isoform_dict))


iter <- 1
while (TRUE) {
  print(paste("Starting iteration",iter))
  start_time <- proc.time()
  log.phi.vi <- phi.update(gamma.vi, tau.vi, log.f.vi, num_pixel, num_celltype, 
                           num_gene, spatial, isoform_dict)
  
  gamma.vi <- gamma.update(alpha.vi, log.phi.vi, num_pixel, num_celltype, spatial)
  
  # check update of tau
  tau.vi <- tau.update(eta*prior_constant_data, log.phi.vi, num_celltype, num_gene,
                       spatial, isoform_dict)
  
  log.f.vi <- f.update(log.phi.vi, num_celltype, num_gene, spatial, isoform_dict)
  
  rownames(log.f.vi) = cell_types
  colnames(log.f.vi) = isoform_vec
  alpha.vi <- alpha.update(alpha.vi, gamma.vi, num_pixel, num_celltype)
  print((proc.time() - start_time)[3])
  iter <- iter + 1
  
  if (iter%%20 != 0) {
    next
  }
  ELBOs <- c(ELBOs, ELBO(eta*prior_constant_data, alpha.vi, gamma.vi, log.phi.vi,  
                         tau.vi, log.f.vi, num_pixel, num_celltype, num_gene, 
                         spatial, isoform_dict))
  print(ELBOs)
  
  print(paste("Alpha length:", length(alpha.vi)))
  print(paste("Gamma dimensions:", dim(gamma.vi)))
  print(paste("Phi dimensions:", dim(log.phi.vi)))
  print(paste("Tau dimensions:",dim(tau.vi)))
  print(paste("f dimensions:",dim(log.f.vi)))
  
  if (ELBOs[length(ELBOs)] - ELBOs[length(ELBOs) - 1] <= 0.00001) {
    break
  }
}

phi.vi = exp(log.phi.vi)
f.vi = exp(log.f.vi)

rownames(gamma.vi) = rownames(spatial)
colnames(gamma.vi) = rownames(reference)
dimnames(phi.vi) <- list( 
  Pixel = rownames(spatial),  # Label for the first dimension (pixels)
  Gene = colnames(spatial),    # Label for the second dimension (genes)
  CellType = rownames(reference)               # Label for the third dimension (cell types)
)
rownames(tau.vi) = rownames(reference)
colnames(tau.vi) = colnames(reference)

output = list(ELBO = ELBOs, eta = eta, alpha = alpha.vi, gamma = gamma.vi,
              phi = phi.vi, tau = tau.vi, f=f.vi, total_iter = iter)

saveRDS(output, file=output_directory)



















