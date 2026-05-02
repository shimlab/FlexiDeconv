#!/bin/bash
#SBATCH --job-name=rctd_v2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data/Joglekar_Reference/logs/rctd_v2.out
#SBATCH --error=/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data/Joglekar_Reference/logs/rctd_v2.err
 
mkdir -p /data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data/Joglekar_Reference/logs
 
module load GCCcore/12.3.0
module load R/4.3.2
 
Rscript /data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Scripts/run_rctd_phase_c.R
 