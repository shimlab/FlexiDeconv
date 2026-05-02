## ============================================================
## Phase C: Run RCTD on V2 Visium HD 8µm Bins
## Assumes Phase A (reference) already completed and saved
## Rebuilds spatial query from Space Ranger output
## ============================================================
library(data.table)
library(Matrix)
library(spacexr)
library(Seurat)
library(nanoparquet)
REF_DIR <- "/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data/Joglekar_Reference"
BASE    <- "/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data"
SR_DIR  <- file.path(BASE, "Visium_Data/spaceranger")
v2_8um_dir <- file.path(SR_DIR, "binned_outputs/square_008um")
## ── Load saved reference from Phase A ──
cat("Loading reference ...\n")
reference <- readRDS(file.path(REF_DIR, "spacexr_reference_no_CB.rds"))
cat("  Done\n\n")
## ── Rebuild spatial query ──
cat("Building spatial query ...\n")
counts_v2 <- Read10X(file.path(v2_8um_dir, "filtered_feature_bc_matrix"))
pos <- as.data.frame(nanoparquet::read_parquet(
  file.path(v2_8um_dir, "spatial", "tissue_positions.parquet")
))
pos <- pos[pos$in_tissue == 1, ]
shared_bc <- intersect(colnames(counts_v2), pos$barcode)
counts_v2 <- counts_v2[, shared_bc, drop = FALSE]
pos <- pos[match(shared_bc, pos$barcode), ]
coords <- data.frame(
  x = pos$pxl_col_in_fullres,
  y = pos$pxl_row_in_fullres,
  row.names = pos$barcode
)
spatial_query <- SpatialRNA(
  coords = coords,
  counts = counts_v2,
  nUMI   = colSums(counts_v2)
)
cat("  Bins:", ncol(counts_v2), "| Genes:", nrow(counts_v2), "\n")
cat("  Median UMI/bin:", median(colSums(counts_v2)), "\n\n")
## ── Run RCTD ──
cat("Creating RCTD object ...\n")
myRCTD <- create.RCTD(
  spatialRNA = spatial_query,
  reference  = reference,
  max_cores  = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
)
cat("Running RCTD in doublet mode ...\n")
myRCTD <- run.RCTD(
  myRCTD,
  doublet_mode = "doublet"
)
rctd_path <- file.path(REF_DIR, "RCTD_results_no_CB.rds")
saveRDS(myRCTD, rctd_path)
cat("  Saved to:", rctd_path, "\n\n")
## ── Extract results ──
cat("Extracting results ...\n")
results_df <- myRCTD@results$results_df
results_df$barcode <- rownames(results_df)
results_df <- merge(results_df, coords, by.x = "barcode", by.y = "row.names")
cat("Classification breakdown:\n")
print(table(results_df$spot_class))
n_total   <- nrow(results_df)
n_singlet <- sum(results_df$spot_class == "singlet")
n_doublet <- sum(results_df$spot_class %in% c("doublet_certain", "doublet_uncertain"))
n_reject  <- sum(results_df$spot_class == "reject")
cat("\nSummary:\n")
cat("  Total bins:", n_total, "\n")
cat("  Singlets:", n_singlet, "(", round(100 * n_singlet / n_total, 1), "%)\n")
cat("  Doublets:", n_doublet, "(", round(100 * n_doublet / n_total, 1), "%)\n")
cat("  Rejected:", n_reject, "(", round(100 * n_reject / n_total, 1), "%)\n")
singlets <- results_df[results_df$spot_class == "singlet", ]
cat("\nSinglet cell types:\n")
print(table(singlets$first_type))
output_path <- file.path(REF_DIR, "RCTD_annotations_no_CB.csv")
write.csv(results_df, output_path, row.names = FALSE)
cat("\nSaved annotations to:", output_path, "\n")

## ── Extract and save cell-type weights ──
## weights is a matrix (pixels × cell types) with estimated proportions
## Useful for downstream soft-assignment analyses, visualisation,
## and cell-type-constrained permutations in Spl-IsoFind
cat("Extracting cell-type weights ...\n")
weights <- myRCTD@results$weights
weights_df <- as.data.frame(as.matrix(weights))
weights_df$barcode <- rownames(weights_df)
weights_df <- merge(weights_df, coords, by.x = "barcode", by.y = "row.names")
weights_path <- file.path(REF_DIR, "RCTD_weights_no_CB.csv")
write.csv(weights_df, weights_path, row.names = FALSE)
cat("  Saved weights to:", weights_path, "\n")
cat("  Dimensions:", nrow(weights_df), "bins x",
    ncol(weights_df) - 3, "cell types\n")  # subtract barcode, x, y

cat("\n=== Done ===\n")