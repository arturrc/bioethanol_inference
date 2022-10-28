#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J uir18 # Job name
#SBATCH -n 4 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000 # Memory in MB
#SBATCH -o tmp/uir18_varcall_%A.out
#SBATCH -e tmp/uir18_varcall_%A.err

. ./vars.sh

# Load modules
module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/apps/gatk-4.1.3.0

echo [$(date)] Job started! UIR_2018 varcall at $chrom. \($SLURM_JOB_ID\) >> batch.log

# Call variants
java -Xmx4g -XX:ParallelGCThreads=4 -jar $GATK4_HOME/gatk-package-4.1.3.0-local.jar \
    Mutect2 \
    -R data/GCF_000146045.2_R64_genomic.fna \
    -I ${scratch}/UIR_2018_1_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_3_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_7_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_9_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_17_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_21_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_27_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_37_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_47_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_57_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_67_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_97_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_109_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_127_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_157_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_167_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_187_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_217_aligned_reads.final.bam \
    -I ${scratch}/UIR_2018_227_aligned_reads.final.bam \
    -O ${scratch}/UIR_2018.${chrom}.vcf.gz \
    -L ${chrom}

echo [$(date)] Job completed! UIR_2018 varcall at $chrom. \($SLURM_JOB_ID\) >> batch.log
