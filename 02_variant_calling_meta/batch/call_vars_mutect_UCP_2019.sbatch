#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J ucp19 # Job name
#SBATCH -n 4 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000 # Memory in MB
#SBATCH -o tmp/ucp19_varcall_%A.out
#SBATCH -e tmp/ucp19_varcall_%A.err

. ./vars.sh

# Load modules
module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/apps/gatk-4.1.3.0

echo [$(date)] Job started! UCP_2019 varcall at $chrom. \($SLURM_JOB_ID\) >> batch.log

# Call variants
java -Xmx4g -XX:ParallelGCThreads=4 -jar $GATK4_HOME/gatk-package-4.1.3.0-local.jar \
    Mutect2 \
    -R data/GCF_000146045.2_R64_genomic.fna \
    -I ${scratch}/UCP_2019_1_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_5_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_13_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_28_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_40_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_54_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_67_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_81_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_95_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_109_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_123_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_137_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_151_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_172_aligned_reads.final.bam \
    -I ${scratch}/UCP_2019_184_aligned_reads.final.bam \
    -O ${scratch}/UCP_2019.${chrom}.vcf.gz \
    -L ${chrom}

echo [$(date)] Job completed! UCP_2019 varcall at $chrom. \($SLURM_JOB_ID\) >> batch.log
