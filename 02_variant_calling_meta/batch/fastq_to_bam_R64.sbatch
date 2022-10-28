#! /bin/bash -e
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J bam # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:30 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 10000 # Memory in MB
#SBATCH -o tmp/gvcf_%A.out
#SBATCH -e tmp/gvcf_%A.err

. ./vars.sh

# Load modules
module load NGmerge/0.3-fasrc01
module load bwa/0.7.15-fasrc02
module load jdk/1.8.0_45-fasrc01
PICARD_HOME=/n/sw/fasrcsw/apps/Core/picard/2.9.0-fasrc01
module load intel/2017c impi/2017.4.239 SAMtools/1.9
module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/apps/gatk-4.1.3.0
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 BCFtools/1.9

echo [$(date)] Job started! Aligning $label. \($SLURM_JOB_ID\) >> batch.log

# Merge reads
NGmerge -1 $fastq1 -2 $fastq2 -a -v -z -o ${scratch}/meta_${label}.trimmed

# Align to s288c genome
cd data
bwa mem -M -t 1 \
    -R "@RG\tID:i${label}\tSM:${label}\tPL:ILLUMINA" \
    GCF_000146045.2_R64_genomic ${scratch}/meta_${label}.trimmed_1.fastq.gz ${scratch}/meta_${label}.trimmed_2.fastq.gz > \
    ${scratch}/${label}_aligned_reads.sam
cd ..

echo [$(date)] Job progress! Alignment on $label done. \($SLURM_JOB_ID\) >> batch.log

# Sort and convert sam to bam file
java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam \
    I=${scratch}/${label}_aligned_reads.sam \
    O=${scratch}/${label}_aligned_reads.sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

# Mark duplicates
java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar MarkDuplicates \
    I=${scratch}/${label}_aligned_reads.sorted.bam \
    O=${scratch}/${label}_aligned_reads.dedup.bam \
    METRICS_FILE=tmp/${label}_aligned_reads.dedup_metrics.txt \
    REMOVE_DUPLICATES=false \
    TAGGING_POLICY=All

# Resort and reindex
java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam \
    I=${scratch}/${label}_aligned_reads.dedup.bam \
    O=${scratch}/${label}_aligned_reads.final.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

# Validate bam file
java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar ValidateSamFile \
    I=${scratch}/${label}_aligned_reads.final.bam \
    O=out/${label}_aligned_reads.validate.txt \
    MODE=SUMMARY

echo [$(date)] Job completed! Alignment of $label. \($SLURM_JOB_ID\) >> batch.log
