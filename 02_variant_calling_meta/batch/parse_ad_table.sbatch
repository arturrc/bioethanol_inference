#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J parse # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:30 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 1000 # Memory in MB
#SBATCH -o tmp/parse_%A.out
#SBATCH -e tmp/parse_%A.err

. ./vars.sh

# Load modules
module load R/4.1.0-fasrc01

echo [$(date)] Job started! Parsing table of $expt. \($SLURM_JOB_ID\) >> batch.log

# Parse table
input=${scratch}/${expt}.table
out_counts=out/${expt}_counts.txt
out_depths=out/${expt}_depths.txt
Rscript parse_ad_table.r $input $out_counts $out_depths

echo [$(date)] Job completed! Parsing table of $expt. \($SLURM_JOB_ID\) >> batch.log
