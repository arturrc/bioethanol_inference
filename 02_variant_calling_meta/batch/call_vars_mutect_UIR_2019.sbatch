#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J uir19 # Job name
#SBATCH -n 4 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000 # Memory in MB
#SBATCH -o tmp/uir19_varcall_%A.out
#SBATCH -e tmp/uir19_varcall_%A.err

. ./vars.sh

# Load modules
module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/apps/gatk-4.1.3.0

echo [$(date)] Job started! UIR_2019 varcall at $chrom. \($SLURM_JOB_ID\) >> batch.log

# Call variants
java -Xmx4g -XX:ParallelGCThreads=4 -jar $GATK4_HOME/gatk-package-4.1.3.0-local.jar \
    Mutect2 \
    -R data/GCF_000146045.2_R64_genomic.fna \
    -I ${scratch}/UIR_2019_1_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_8_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_13_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_27_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_38_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_51_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_65_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_77_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_92_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_107_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_119_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_134_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_150_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_163_aligned_reads.final.bam \
    -I ${scratch}/UIR_2019_174_aligned_reads.final.bam \
    -O ${scratch}/UIR_2019.${chrom}.vcf.gz \
    -L ${chrom}

echo [$(date)] Job completed! UIR_2019 varcall at $chrom. \($SLURM_JOB_ID\) >> batch.log
