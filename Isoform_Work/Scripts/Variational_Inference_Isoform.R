library(bayess)



# Read spatial data and reference
#spatial <- readRDS("../../Data/Isoform/Simulation/spatial_data.rds")
spatial <- readRDS("../../Data/Isoform/LebrigandIsoform/proc_data/iso_LR_spatial.rds")
#reference <- readRDS("../../Data/Isoform/Simulation/reference.rds")
reference <- readRDS("../../Data/Isoform/LebrigandIsoform/proc_data/iso_LR_ref.rds")
#isoform_dict <- readRDS("../../Data/Isoform/Simulation/isoform_dict.rds")
isoform_dict <- readRDS("../../Data/Isoform/LebrigandIsoform/proc_data/iso_dict.rds")

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
#prior_constant_data = readRDS(prior_constant)
prior_constant_data = readRDS("../../Data/Isoform/LebrigandIsoform/proc_data/iso_LR_pw.rds")
medium_molecule_count = median(rowSums(spatial))
prior_constant_data = prior_constant_data * medium_molecule_count * num_pixel/num_celltype


# Define Output directory
#output_directory <- args[4]



################################################################################
source("/data/projects/punim0614/yichenjiang/LDAproject/Scripts/VI_Rscripts_Isoform.R")


alpha.vi <- rdirichlet(n = 1, par = rep(5, num_celltype))[1,] * 50
gamma.vi <- num_celltype * rdirichlet(n = num_pixel, par = rep(30, num_celltype))

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
log.phi.vi = log(phi.vi)
eta <- reference


ELBOs <- c(ELBO(eta*prior_constant_data, gamma.vi, tau.vi, alpha.vi, phi.vi, num_celltype, 
                f.vi, num_pixel, num_gene, spatial, isoform_dict))


iter <- 1
while (TRUE) {
  print(paste("Starting iteration",iter))
  #print(paste("Alpha length:", length(alpha.vi)))
  #print(paste("Gamma dimensions:", dim(gamma)))
  #print(paste("Phi dimensions:", dim(phi)))
  #print(paste("Tau dimensions:",dim(tau)))
  #print(paste("Eta dimensions:",dim(eta)))
  start_time <- proc.time()
  phi.vi <- phi.update(num_pixel, num_celltype, num_gene, gamma.vi, tau.vi, 
                    isoform_dict, f.vi, num_transcript)
  
  gamma.vi <- gamma.update(num_pixel, num_celltype, alpha.vi,phi.vi, 
                        spatial)
  
  # check update of tau
  tau.vi <- tau.update(num_celltype, num_pixel, eta*prior_constant_data, phi.vi, spatial,
                    num_transcript, isoform_dict, num_gene)
  
  f.vi <- f.update(num_celltype, isoform_dict, num_gene, phi.vi, spatial)
  
  
  
  rownames(f.vi) = cell_types
  colnames(f.vi) = isoform_vec
  alpha.vi <- alpha.update(num_pixel, num_celltype, alpha.vi, gamma.vi)
  print((proc.time() - start_time)[3])
  iter <- iter + 1
  
  if (iter%%50 != 0) {
    next
  }
  ELBOs <- c(ELBOs, ELBO(eta*prior_constant_data, gamma.vi, tau.vi, alpha.vi, phi.vi, num_celltype, 
                         f.vi, num_pixel, num_gene, spatial, isoform_dict))
  print(ELBOs)
  if (ELBOs[length(ELBOs)] - ELBOs[length(ELBOs) - 1] <= 0.00001) {
    break
  }
}


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
#target.dir <- paste(DATA_DIR, "/VI_result/result_", PRIOR_CONST, "_", file_trail
#                    , ".rds", sep="")
saveRDS(output, file="../../Data/Isoform/Simulation/results/initialRes/initialOutput.rds")









