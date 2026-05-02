library(STdeconvolve)
library(Matrix)
source("/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Scripts/functions(new).R")

args <- commandArgs(trailingOnly = TRUE)


spatial <- args[1]
reference <- args[2]
print(spatial)
print(reference)

st_data = readRDS(spatial)
pixel_num = dim(st_data)[1]
gene_num = dim(st_data)[2]
ref_data = readRDS(reference)
cell_types = rownames(ref_data)

medium_molecule_count = median(rowSums(st_data))



prior_constant <- args[3]
if (is_number(prior_constant)) {
  PRIOR_CONST = matrix(as.numeric(prior_constant), nrow = dim(ref_data)[1],
                       ncol = dim(ref_data)[2], byrow = FALSE)
  input_pw = prior_constant
} else {
  # a directory
  prior_constant_data = readRDS(prior_constant) * as.numeric(args[5])
  input_pw = prior_constant_data
  PRIOR_CONST = matrix(prior_constant_data, nrow = dim(ref_data)[1],
                       ncol = dim(ref_data)[2], byrow = FALSE)
  PRIOR_CONST = PRIOR_CONST * medium_molecule_count * pixel_num/length(cell_types)
}
print(PRIOR_CONST[,1])

output_directory <- args[4]

#spatial="/data/gpfs/projects/punim0614/yichen/Data/STdecon_sim/MERFISHMPOA/spatial_data.rds"
#reference="/data/gpfs/projects/punim0614/yichen/Data/STdecon_sim/MERFISHMPOA/reference.rds"
#prior_constant=500
#output_directory="/data/gpfs/projects/punim0614/yichen/Data/STdecon_sim/MERFISHMPOA/results/pc500.rds"


print(output_directory)

# This file requires Spatial Transcriptomics dataset and Reference Dataset

# Spatial Transcriptomics dataset should be gene x spots
# Reference is cell type x gene
#DATA_DIR = "/data/projects/punim0614/yichen/Data"
#version = paste("/svg_data_", file_trail, ".rds", sep="")



#ref_data = readRDS(paste(DATA_DIR, "/svg_data_ref.rds", sep=""))


print(paste("Dimension of Spatial Data: ", dim(st_data)))
print(paste("Dimension of Reference Data: ", dim(ref_data)))


###############################################################################

# Randomly sample 196 spots to visualize
# gt_file_path = paste(DATA_DIR, "/cell_type_ground_truth_", file_trail, ".rds", sep="")
# st_ground_truth = readRDS(gt_file_path)
# num_spots = 196
#
# spot_idx = sort(sample(1:length(st_ground_truth$cell_type_one), size=num_spots))
# SIDE.LENGTH <- sqrt(num_spots)
# center.x <- rep(seq(1,SIDE.LENGTH)*50,SIDE.LENGTH)
# center.y <- rep(seq(1,SIDE.LENGTH)*50,each=SIDE.LENGTH)
# positions <- data.frame(x = center.x, y = center.y)
#
#
# cell_type_gt <- matrix(0, nrow = num_spots, ncol = length(cell_types))
# colnames(cell_type_gt) = cell_types
# for (i in 1:num_spots) {
#   ct1 = st_ground_truth$cell_type_one[spot_idx[i]]
#   ct2 = st_ground_truth$cell_type_two[spot_idx[i]]
#   prop1 = st_ground_truth$cell_type_prop1[spot_idx[i]]
#   prop2 = st_ground_truth$cell_type_prop2[spot_idx[i]]
#   cell_type_gt[i,ct1] = prop1
#   cell_type_gt[i,ct2] = prop2
# }
# vizAllTopics(cell_type_gt, positions, r=20)

###############################################################################
source("/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Scripts/VI_Rscripts.R")



alpha.vi <- rdirichlet(n = 1, par = rep(5, length(cell_types)))[1,] * 50
#alpha.vi <- rep(50/length(cell_types), length(cell_types))
print(paste("Alpha length:", length(alpha.vi)))

# gamma is pixel x cell type
gamma <- length(cell_types) * rdirichlet(n = pixel_num, par = rep(30, length(cell_types)))
#gamma <- matrix(1, nrow = pixel_num, ncol = length(cell_types))

print(paste("Gamma dimensions:", dim(gamma)))


# phi is pixel x gene x cell type
phi <- array(0, dim = c(pixel_num, gene_num, length(cell_types)))
for (i in 1:pixel_num) {
  phi[i,,] = rdirichlet(n = gene_num, par = rep(30, length(cell_types)))
}
#phi <- array(1/length(cell_types), dim = c(pixel_num, gene_num, length(cell_types)))

print(paste("Phi dimensions:", dim(phi)))


# tau is cell type x gene
tau <- length(cell_types) * rdirichlet(n = length(cell_types), par = rep(10, gene_num))
#tau <- matrix(1, nrow = length(cell_types), ncol = gene_num)

print(paste("Tau dimensions:",dim(tau)))

eta <- ref_data
print(paste("Eta dimensions:",dim(eta)))


ELBOs <- c(ELBO(eta*PRIOR_CONST, gamma, tau, alpha.vi, phi, length(cell_types),
                pixel_num, gene_num, st_data))
iter <- 1
while (TRUE) {
  print(paste("Starting iteration",iter))
  #print(paste("Alpha length:", length(alpha.vi)))
  #print(paste("Gamma dimensions:", dim(gamma)))
  #print(paste("Phi dimensions:", dim(phi)))
  #print(paste("Tau dimensions:",dim(tau)))
  #print(paste("Eta dimensions:",dim(eta)))
  start_time <- proc.time()
  phi <- phi.update.old(pixel_num, length(cell_types), gene_num, gamma, tau)
  gamma <- gamma.update(pixel_num, gene_num, length(cell_types), alpha.vi, phi,
                        st_data)
  tau <- tau.update(length(cell_types), gene_num, pixel_num,
                    eta*PRIOR_CONST, phi, st_data)
  alpha.vi <- alpha.update(pixel_num, length(cell_types), alpha.vi, gamma)
  print((proc.time() - start_time)[3])
  iter <- iter + 1

  if (iter%%50 != 0) {
    next
  }
  ELBOs <- c(ELBOs, ELBO(eta*PRIOR_CONST, gamma, tau, alpha.vi, phi, length(cell_types),
                         pixel_num, gene_num, st_data))
  print(ELBOs)
  if (ELBOs[length(ELBOs)] - ELBOs[length(ELBOs) - 1] <= 0.0001) {
    break
  }
}

rownames(gamma) = rownames(st_data)
colnames(gamma) = rownames(ref_data)
dimnames(phi) <- list(
  Pixel = rownames(st_data),  # Label for the first dimension (pixels)
  Gene = colnames(st_data),    # Label for the second dimension (genes)
  CellType = rownames(ref_data)               # Label for the third dimension (cell types)
)
rownames(tau) = rownames(ref_data)
colnames(tau) = colnames(ref_data)

output = list(ELBO = ELBOs, eta = eta, alpha = alpha.vi, gamma = gamma,
              phi = phi, tau = tau, total_iter = iter, prior_weight = input_pw)
#target.dir <- paste(DATA_DIR, "/VI_result/result_", PRIOR_CONST, "_", file_trail
#                    , ".rds", sep="")
saveRDS(output, file=output_directory)
print("Saved output")
