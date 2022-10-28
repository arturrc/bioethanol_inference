# METAGENOMIC VARIANT CALLING SCRIPTS -------------------------------------------------------------------------------
Scripts for calling variants from metagenomic reads.
This pipeline is designed to run on a machine running CentOS 7 with Slurm job management.
This means that significant adaptation in terms of calling dependencies and running job
scripts will be necessary.

# Dependencies
# Specific versions are sufficient, but not necessarily required.
	- bwa 0.7.15
	- SAMtools 1.9
	- GATK 4.0.2.1
	- JDK 1.8.0_45
	- NGmerge 0.3
	- BCFtools 1.9
	- R 4.1.0

# Structure
	- The root ./ (where this README file is) contains all working scripts.
	- All scripts should be run from the root ./ .
	- ./runall.sh guides through the whole pipeline.
	- ./batch/ contains Slurm job scripts.
	- ./data/ contains the s288c reference genome from SGD (https://www.yeastgenome.org/),
	and should contain read files or links to read files (stored in NCBI Bioproject PRJNA865262).
