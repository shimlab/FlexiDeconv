#!/bin/bash
#SBATCH --job-name=Lebrigand_isoform_logsumtrick
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=47:00:00               # Time limit hrs:min:sec
#SBATCH --output=/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Outputs/LebIso_LogSum_Output.out     # Standard output
#SBATCH --error=/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Outputs/LebIso_LogSum_Output.err      # Standard error log file

START_TIME_=$(date +%s)
st=$(date)

module load R/4.4.0


vi_file="/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Scripts/Variational_Inference_Isoform.R"

Rscript "$vi_file"


# Adding the startime
echo
echo "Running on: $(hostname)"
echo "Start time: $st"

# Adding the endtime
END_TIME_=$(date +%s)
et=$(date)
echo "End time: $et"

runtime=$((END_TIME_ - START_TIME_))

# Append the runtime to the output file
echo "Runtime: $runtime seconds"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
