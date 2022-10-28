# LINEAGE INFERENCE PIPELINE
# Throughout the code: 
# - expt refers to site-year,
# - cluster refers to inferred lineages of descent,
# - geno refers to clonal isolate sequencing data,
# - meta refers to metagenomic sequencing data,
# - diagnostic vars refers to lineage-specific variants
#   used to infer lineage frequency.

# 0) Set up output folder
mkdir -p out

# 1) Assign clusters based on clone phylogenies
# Reroot Tree 1 based on where bioethanol strains are rooted in Tree 2.
Rscript reroot_tree.R
# Assign clusters based on rerooted tree
Rscript assign_clusters.R

# 2) Identify cluster diagnostic variants from geno data
Rscript find_diagnostic_vars.R
# Count the number of diagnostic variants overlapping between geno and meta
Rscript count_overlap_diag_vars.R

# 3) Estimate lineage freqs independently for each cluster
# Useful to QC the joint inference done later
Rscript estimate_lineage_freqs_independently.R

# 4) Estimate lineage freqs jointly
# This is done on the cluster. Use ./runall_cluster.sh to carry
# out this portion of the pipeline.

# 5) Plot Mullers and trees
Rscript choose_lineage_colors.R
Rscript plot_labeled_tree.R
Rscript plot_tree_ploidy.R
Rscript plot_circ_tree.R
Rscript plot_mullers.R

# 6) Calculate frequency in population (as opposed to frequency in the metagenome)
# and plot Mullers
Rscript calc_popfreqs.R
Rscript plot_mullers_popfreq.R

# EXTRA ==================================================
# Plot allele frequencies and coverage for each clone along the genome
Rscript plot_allele_freq_and_coverage.R
Rscript plot_coverage_rank.R

# Plot subset of raw metagenomic data
Rscript plot_meta_lines.R

# Plot raw metagenomic data for diagnostic variants for validation
Rscript plot_diagnostic_vars_trajs.R

# Plot distributon of maximum inferred frequency per lineage
Rscript plot_dist_max_inferred_freq.R

# Plot genotype posterior probability for different counts, depths, and error rates
Rscript plot_post_geno_prob_examples.R

# Interactive script that counts mutations and overlaps over different datasets
# characterize_diversity.R