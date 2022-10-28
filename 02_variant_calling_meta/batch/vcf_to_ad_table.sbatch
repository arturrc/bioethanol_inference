#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J genotype # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:30 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 5000 # Memory in MB
#SBATCH -o tmp/genotype_%A.out
#SBATCH -e tmp/genotype_%A.err

. ./vars.sh

# Load modules
module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/apps/gatk-4.1.3.0

echo [$(date)] Job started! Genotyping for $chrom. \($SLURM_JOB_ID\) >> batch.log

# Make database
java -Xmx4g -XX:ParallelGCThreads=1 -jar $GATK4_HOME/gatk-package-4.1.3.0-local.jar \
    VariantsToTable \
    -V ${vcf} \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF AD \
    -O ${output}

echo [$(date)] Job completed! Genotyped at $chrom. \($SLURM_JOB_ID\) >> batch.log

