# 1) READ TRIMMING
module load NGmerge/0.2-fasrc01

for i in 10 23 24; do j=$(expr ${i} + 113); NGmerge -1 AKG${I}_S${J}_R1_001.1.fastq -2 AKG${I}_S${J}_R2_001.2.fastq -a -v -o AKG${I}_S${J}.trimmed; done

# 2) ALIGNMENT AGAINST REFERENCE GENOME
module load bwa/0.7.15-fasrc02

# indexing the ref. genome with BWA needs to be done once
# this will create 5 additional files
bwa index -p GCF_000146045.2_R64_genomic GCF_000146045.2_R64_genomic.fna

# mem is the BWA function for alignment
for i in 10 23 24; do j=$(expr ${i} + 113); bwa mem -M -t 1 -R "@RG\tID:HTLMC.1\tSM:AKG${i}_S${j}\tPL:ILLUMINA" GCF_000146045.2_R64_genomic AKG${i}_S${j}.trimmed_1.fastq AKG${i}_S${j}.trimmed_2.fastq > AKG${i}_S${j}.sam 2> AKG${i}_S${j}.bwa.log; done

# 3) CONVERSION FROM SAM TO BAM
# can be done with Samtools or Picard
# here Picard is used and data are also sorted and indexed
# Picard runs on java

module load jdk/1.8.0_45-fasrc01
PICARD_HOME=/n/sw/fasrcsw/apps/Core/picard/2.9.0-fasrc01

for i in 10 23 24; do j=$(expr ${i} + 113); java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam \
	I=AKG${i}_S${j}.sam \
	O=AKG${i}_S${j}.sorted.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true
done

# 4) MARKING DUPLICATES
# also performed with Picard

for i in 10 23 24; do j=$(expr ${i} + 113); java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar MarkDuplicates \
	I=AKG${i}_S${j}.sorted.bam \
	O=AKG${i}_S${j}.dedup.bam \
	METRICS_FILE=AKG${i}_S${j}.dedup_metrics.txt \
	REMOVE_DUPLICATES=false \
	TAGGING_POLICY=All
done 2> AKG10_23and24_dedup.log

# 5) RESORTING AND REINDEXING
# SortSam (Picard) needs to be run again

for i in 10 23 24; do j=$(expr ${i} + 113); java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam \
	I=AKG${i}_S${j}.dedup.bam \
	O=AKG${i}_S${j}.final.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true
done 2> final_sorting.log

# 6) VALIDATING THE BAM FILES
# also performed with Picard

for i in 10 23 24; do j=$(expr ${i} + 113); java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar ValidateSamFile \
	I=AKG${i}_S${j}.final.bam \
	O=AKG${i}_S${j}.validate.txt \
	MODE=SUMMARY
done

# 7) IF "no errors found", ALL PREVIOUS SAM AND BAM FILES CAN NOW
# BE DELETED (IN PRINCIPLE IT IS ALWAYS POSSIBLE TO REGENERATE
# EVEN THE ORIGINAL FASTQ FILES FROM THE FINAL BAM FILE)

# 8) HAPLOTYPE CALLING
# we will use GATK4, the most recent version (Gold standard)
# 3 steps are needed: haplotype calling, database creation
# and genotyping (this is 1 step more than with GATK3)
# GATK4 needs to be copied to a known directory

module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/gatk-4.0.8.1

# 8A) generating the g.vcf files for the cohort samples

for i in 10 24; do j=$(expr ${i} + 113); java -Xmx4g -XX:ParallelGCThreads=1 -jar $GATK4_HOME/gatk-package-4.0.8.1-local.jar HaplotypeCaller -R ~/genome_references/GCF_000146045.2_R64_genomic.fna -I /n/desai_lab/users/andygombert/output_Bwa/AKG${i}_S${j}.final.bam -O /n/desai_lab/users/andygombert/output_Bwa/AKG${i}_S${j}.g.vcf --emit-ref-confidence GVCF; done

# options below only used for low sequence coverage (< 10x)
# --min-pruning 1
# --min-dangling-branch-length 1

# gatk_java="-Xmx4g -XX:ParallelGCThreads=1"
# I should request the at least that memory (4G) and this no. of threads+1

# 8B) generating a databse from the g.vcf files
# here all chromosomes are listed; can be chosen at will

~/gatk-4.0.8.1/gatk GenomicsDBImport -R ~/genome_references/GCF_000146045.2_R64_genomic.fna -V /n/desai_lab/users/andygombert/output_Bwa/AKG10_S123.g.vcf -V /n/desai_lab/users/andygombert/output_Bwa/AKG23_S136.g.vcf -V /n/desai_lab/users/andygombert/output_Bwa/AKG24_S137.g.vcf --genomicsdb-workspace-path database_032819 \
-L NC_001133.9 \
-L NC_001134.8 \
-L NC_001135.5 \
-L NC_001136.10 \
-L NC_001137.3 \
-L NC_001138.5 \
-L NC_001139.9 \
-L NC_001140.6 \
-L NC_001141.2 \
-L NC_001142.9 \
-L NC_001143.9 \
-L NC_001144.5 \
-L NC_001145.3 \
-L NC_001146.8 \
-L NC_001147.6 \
-L NC_001148.4 \
-L NC_001224.1

# 8C) genotyping from the database

~/gatk-4.0.8.1/gatk GenotypeGVCFs -R ~/genome_references/GCF_000146045.2_R64_genomic.fna -V gendb://database_032819 -O /n/desai_lab/users/andygombert/output_Bwa/AKG10_23and24.vcf --heterozygosity 0.005 2> /n/desai_lab/users/andygombert/output_Bwa/AKG10_23and24_genotyped.log
