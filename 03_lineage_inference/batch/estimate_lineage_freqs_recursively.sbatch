#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J infer # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:30 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 1000 # Memory in MB
#SBATCH -o tmp/infer_%A.out
#SBATCH -e tmp/infer_%A.err

module load R/3.3.3-fasrc01

# Import variables
. ./vars.sh

echo [$(date)] Job started! Estimating lineage freqs for ${expt}, tp \#${tp_id}. \($SLURM_JOB_ID\) >> batch.log

# Run inference
Rscript estimate_lineage_freqs_recursively.r ${expt} ${tp_id} ${scratch}/freqs_${expt}_${tp_id}.txt ${scratch}/fit_stats_${expt}_${tp_id}.txt

echo [$(date)] Job completed! Estimating lineage freqs for ${expt}, tp \#${tp_id}. \($SLURM_JOB_ID\) >> batch.log
