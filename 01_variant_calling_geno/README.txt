# CLONAL ISOLATE VARIANT CALLING SCRIPTS -------------------------------------------------------------------------------
Scripts for calling variants from clonal isolate sequencing data.
This pipeline is designed to run on a machine running CentOS 7 with Slurm job management.
This means that significant adaptation in terms of calling dependencies and running job
scripts will be necessary.
Input sequencing read files are ultimately stored in NCBI Bioproject PRJNA865262.
The pipeline commands are laid out in ./runall.sh.
./GCF_000146045.2_R64_genomic.fna is the s288c reference genome from SGD (https://www.yeastgenome.org/).

# Dependencies
# Specific versions are sufficient, but not necessarily required.
	- bwa 0.7.15
	- SAMtools 1.9
	- GATK 4.0.2.1
	- JDK 1.8.0_45
	- NGmerge 0.3
	- BCFtools 1.9
	- R 4.1.0
