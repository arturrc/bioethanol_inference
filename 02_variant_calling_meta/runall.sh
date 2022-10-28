. ./vars.sh

# 0) Create directories project
mkdir -p $scratch
mkdir -p $outdir
mkdir -p out
mkdir -p tmp

# 0) Initiate log file
touch batch.log

# 0) Import read data
# In this case, files are linked from some external location.
# Read files are ultimately available in NCBI Bioproject PRJNA865262.
. ./link_data.sh

# 0) Index reference genome for bwa
module load bwa/0.7.15-fasrc02
cd data
bwa index -p GCF_000146045.2_R64_genomic GCF_000146045.2_R64_genomic.fna
cd ..

# 0) Prepare reference genome for GATK
module load intel/2017c impi/2017.4.239 SAMtools/1.9
samtools faidx data/GCF_000146045.2_R64_genomic.fna

module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/apps/gatk-4.1.3.0
java -Xmx4g -XX:ParallelGCThreads=1 -jar $GATK4_HOME/gatk-package-4.1.3.0-local.jar \
    CreateSequenceDictionary \
    -R data/GCF_000146045.2_R64_genomic.fna

# PIPELINE =========================================================================
# 1) Align reads to s288c genome
for i in {1..69}
do
    echo \> ${i}/69
    label=${meta_labels[i]}
    code=${meta_fastqgz_labels[i]}
    fastq1=data/${code}_*_R1_*.fastq.gz
    fastq2=data/${code}_*_R2_*.fastq.gz
    sbatch --export=fastq1=${fastq1},fastq2=${fastq2},label=${label} batch/fastq_to_bam_R64.sbatch
    sleep 0.2
done

# 2) Call variants by site-year, and by chromosome
for expt in ${expts[@]}
do
    for i in {0..16}
    do
        echo \> ${i}/16
        chrom=${chroms[i]}
        sbatch --export=chrom=${chrom} batch/call_vars_mutect_${expt}.sbatch
        sleep 0.2
    done
done

# 3) Extract AD to table
for expt in ${expts[@]}
do
    for i in {0..16}
    do
        echo \> ${i}/16
        chrom=${chroms[i]}
        sbatch --export=expt=${expt},chrom=${chrom} batch/vcf_to_table.sbatch
        sleep 0.2
    done
done

# 4) Merge tables, skip mitochondria
for expt in ${expts[@]}
do
    chrom=${chroms[0]}
    cat ${scratch}/${expt}.${chrom}.table > ${scratch}/${expt}.table
    for i in {1..15}
    do
        echo \> ${expt}: ${i}/16
        chrom=${chroms[i]}
        tail -n +2 ${scratch}/${expt}.${chrom}.table >> ${scratch}/${expt}.table
    done
done

# 5) Parse AD tables
for expt in ${expts[@]}
do
    echo \> ${expt}
    sbatch --export=expt=${expt} batch/parse_ad_table.sbatch
    sleep 0.2
done
