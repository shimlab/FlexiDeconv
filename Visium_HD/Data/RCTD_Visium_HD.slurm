#!/bin/bash
#SBATCH --job-name=rctd_v2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data/Joglekar_Reference/logs/rctd_no_CB.out
#SBATCH --error=/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data/Joglekar_Reference/logs/rctd_no_CB.err
 
mkdir -p /data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data/Joglekar_Reference/logs
 
module load GCC/11.3.0 OpenMPI/4.1.4
module load R/4.4.0
 
Rscript /data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Visium_HD/Data/RCTD_Visium_HD.R
 