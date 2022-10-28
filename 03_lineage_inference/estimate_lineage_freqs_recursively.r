rm(list = ls())
options(stringsAsFactors = F)

# PREAMBLE =============================================
library(reshape2)
library(plyr)
library(matrixStats)
library(R.utils)

source("vars.R")
source("functions.R")

args <- commandArgs(trailingOnly = T)
expt <- args[1]
tp_id <- as.numeric(args[2])
output_freqs <- args[3]
output_fit_stats <- args[4]

# CODE =================================================
cat(c("# ", expt, "\n"), sep = "")

# Get relevant experiment info
all_clones <- c(expt_info[[expt]]$starter_clones, expt_info[[expt]]$picked_clones)
meta_tps <- expt_info[[expt]]$meta_tps

# Import cluster info
clusters <- read.table(sprintf("out/clusters_%s.txt", expt), header = T, sep = "\t", check.names = F)
clusters <- sapply(1:nrow(clusters), function(i) strsplit(clusters$clones[i], ";")[[1]], simplify = F)
cluster_sizes <- sapply(clusters, length)

# Import inference cluster info
cluster_info <- read.table(sprintf("out/cluster_fit_info_%s.txt", expt), sep = "\t", header = T)

# Import diagnostic variants 
diag_vars <- read.table(sprintf("out/diagnostic_vars_%s.txt", expt), header = T, sep = "\t", check.names = F)
diag_vars <- sapply(seq_along(clusters), function(i) {
  if(i %in% diag_vars$cluster){
    out <- diag_vars[diag_vars$cluster == i,]
    out_geno <- do.call(rbind, strsplit(out$geno_log_posts, ";"))
    out_geno <- matrix(as.numeric(out_geno), ncol = ncol(out_geno))
    out_geno <- as.data.frame(out_geno)
    names(out_geno) <- paste0("geno_log_post_", 0:(ncol(out_geno) - 1))
    out <- cbind(out[,1:6], out_geno)
    return(out)
  } else {
    return(NA)
  }
}, simplify = F)

# Read geno count and depth data
geno_counts <- read.table(sprintf("data/geno/geno_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
geno_depths <- read.table(sprintf("data/geno/geno_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
geno_n_vars <- nrow(geno_counts)

# Read meta count and depth data
meta_counts <- read.table(sprintf("data/%s_counts.txt", expt), sep = "\t", header = T, check.names = F)
meta_depths <- read.table(sprintf("data/%s_depths.txt", expt), sep = "\t", header = T, check.names = F)
meta_n_vars <- nrow(meta_counts)

# Find diagnostic variant indices in both genotype and metagenome data
diag_vars_overlap <- geno_meta_overlap(diag_vars, geno_counts, meta_counts)
geno_diag_ids <- diag_vars_overlap$geno_ids
meta_diag_ids <- diag_vars_overlap$meta_ids

# Create list of genotypes for all clusters
genotypes <- sapply(seq_along(clusters), function(i){
  cluster_clones <- clusters[[i]][1]
  ploidy <- ploidies[[cluster_clones[1]]]
  out <- seq(0, 1, 1/ploidy)
  return(out)
}, simplify = F)

# Get boolean flip vectors for overlapping variants
diag_flip <- sapply(seq_along(clusters), function(i){
  x <- geno_diag_ids[[i]]
  if(is.na(list(x))){
    return(NA)
  } else {
    out <- sapply(x, function(y) diag_vars[[i]]$flip[diag_vars[[i]]$diag_var == y])
    return(out)
  }
}, simplify = F)

# Get geno_log_posts for overlapping variants
geno_log_posts <- sapply(seq_along(clusters), function(i){
  x <- geno_diag_ids[[i]]
  if(is.na(list(x))){
    return(NA)
  } else {
    out <- unname(as.matrix(diag_vars[[i]][match(x, diag_vars[[i]]$diag_var),-c(1:6)]))
    return(out)
  }
}, simplify = F)

# Define clusters to be fit: remove NAs
included_clusters <- which(sapply(geno_diag_ids, function(x) !is.na(x[1])))

# Filter included clusters
ids <- sapply(included_clusters, function(x) cluster_info$keep[cluster_info$cluster == x])
included_clusters <- included_clusters[ids]

# Define parent-child relationships among included clusters only
edges <- list_parents(clusters[included_clusters], included_clusters)

# Infer frequencies recursively
fits <- infer_freqs_recursively(meta_counts, meta_depths, meta_diag_ids, diag_flip, geno_log_posts, genotypes,
                                edges, included_clusters, clusters, meta_tps, tp_ids = tp_id,
                                constrOptim_timeout = 960, # can often be increased if inference fails and times out
                                n_its = 10, verbose = T, debug = F)

# Save output
write.table(fits$freqs, output_freqs, sep = "\t",
            quote = F, row.names = F)

write.table(fits$fit_stats, output_fit_stats, sep = "\t",
            quote = F, row.names = F)
