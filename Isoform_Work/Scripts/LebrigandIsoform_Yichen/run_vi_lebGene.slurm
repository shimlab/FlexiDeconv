#!/bin/bash
#SBATCH --job-name=Lebrigand_ISO
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=09:30:00               # Time limit hrs:min:sec
#SBATCH --output=/data/projects/punim0614/yichenjiang/LDAproject/Outputs/lebigrandLR.out     # Standard output
#SBATCH --error=/data/projects/punim0614/yichenjiang/LDAproject/Outputs/lebigrandLR.err      # Standard error log file

START_TIME_=$(date +%s)
st=$(date)

module load R/4.4.0

spatial="/data/gpfs/projects/punim0614/yichenjiang/LDAproject/Data/Isoform/LebrigandIsoform/proc_data/gene_LR_spatial.rds"
reference="/data/gpfs/projects/punim0614/yichenjiang/LDAproject/Data/Isoform/LebrigandIsoform/proc_data/gene_LR_ref.rds"
prior_constant="/data/gpfs/projects/punim0614/yichenjiang/LDAproject/Data/Isoform/LebrigandIsoform/proc_data/gene_LR_pw.rds"
output_dir="/data/gpfs/projects/punim0614/yichenjiang/LDAproject/Data/Isoform/LebrigandIsoform/results/gene_LR_t5.rds"

vi_file="/data/gpfs/projects/punim0614/yichenjiang/LDAproject/Scripts/Variational_Inference.R"

Rscript "$vi_file" "$spatial" "$reference" "$prior_constant" "$output_dir"


# Adding the startime
echo
echo "Start time: $st"

# Adding the endtime
END_TIME_=$(date +%s)
et=$(date)
echo "End time: $et"

runtime=$((END_TIME_ - START_TIME_))

# Append the runtime to the output file
echo "Runtime: $runtime seconds"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
