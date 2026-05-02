#!/bin/bash
#SBATCH --job-name=Lebrigand_isoform_logsumtrick
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=47:00:00               # Time limit hrs:min:sec
#SBATCH --output=/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Outputs/report_results/pw_1.0/LebIso_LogSpace_pw_1.0_assistGamma.out     # Standard output
#SBATCH --error=/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Outputs/report_results/pw_1.0/LebIso_LogSpace_pw_1.0_assistGamma.err      # Standard error log file

START_TIME_=$(date +%s)
st=$(date)

module load R/4.4.0


spatial="/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Data/Isoform/LebrigandIsoform/proc_data_pipeline/iso_LR_spatial.rds"
reference="/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Data/Isoform/LebrigandIsoform/proc_data_pipeline/iso_LR_ref.rds"
isoform_dict="/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Data/Isoform/LebrigandIsoform/proc_data_pipeline/iso_LR_dict.rds"
prior_constant="/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Data/Isoform/LebrigandIsoform/proc_data_pipeline/iso_LR_pw.rds"
output_dir="/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Outputs/report_results/pw_1.0/Output_LogSpace_pw_1.0_assistGamma.rds"
prior_multiplier=10
seed=5

vi_file="/data/gpfs/projects/punim0614/ashtonl/FlexiDeconv/Isoform_Work/Scripts/Variational_Inference_Isoform.R"

Rscript "$vi_file" "$spatial" "$reference" "$isoform_dict" "$prior_constant" "$prior_multiplier" "$output_dir" "$seed"


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
