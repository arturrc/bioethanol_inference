# LINEAGE INFERENCE SCRIPTS ---------------------------------------------------------------------------------------
Scripts for lineage inference and figure generation.

# Dependencies
# Specific versions are sufficient, but not necessarily required.
	- R 4.2.1
	- R packages
		- ape 5.6-2
		- ggmuller 0.5.4
		- ggplot2 3.3.6
		- matrixStats 0.62.0
		- plyr 1.8.7
		- R.utils 2.12.0
		- reshape2 1.4.4

# Structure
	- The root ./ (where this README file is) contains all working scripts.
	- All scripts should be run from the root ./ .
	- ./runall.sh guides through the whole pipeline.
	- The pipeline can be run mostly locally, except for the lineage inference itself, for which 
	  the code is written for a machine running CentOS 7 with Slurm job management (see ./runall.sh).
	- ./batch/ contains Slurm job scripts.
	- ./data/ contains all input data necessary (see description below).

# Data files
	- ./data/geno/ : Files with alternate allele depth and counts from clonal sequencing data split by site-year.
	- ./data/meta/ : Files with alternate allele depth and counts from whole-genome metagenomic sequencing data 
	  split by site-year.
	- ./data/tree1.ml.tree : Newick format of Tree 1. See Methods.
	- ./data/info_expts.txt : Table with information on all site-years analyzed.
	- ./data/info_metagenomics.txt : Table with information on all metagenomics timepoints analyzed.
	- ./data/info_picked_clones.txt : Table with information on all picked fermentation clones.
	- ./data/info_starter_clones.txt : Table with information on all picked starter clones.
	- ./data/info_starter_strains.txt : Table with information on all starter strains used in each site-year.
	- ./data/1011_varids.txt : Table with chromosome, reference allele, and alternate allele from the YGP collection 
	  (Peter et al., Nature v. 556, pp. 339–344, 2018; http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz)
	- ./data/saccharomyces_cerevisiae_R64-3-1_20210421.gff : s288c reference genome annotation from SGD
  	  (https://www.yeastgenome.org/)
