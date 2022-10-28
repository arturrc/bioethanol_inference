rm(list = ls())
options(stringsAsFactors = F)

# PREAMBLE =============================================
library(reshape2)
library(plyr)
library(ggplot2)
library(matrixStats)

source("vars.R")
source("functions.R")

# CODE =================================================
# Set lineage filter
max_zero_depths <- 0
min_n_diag_vars <- 3
min_median_depth <- 10

# Cycle throuogh experiments
for(expt in expts){
  cat(c("# ", expt, "\n"), sep = "")
  all_clones <- c(expt_info[[expt]]$starter_clones, expt_info[[expt]]$picked_clones)
  meta_tps <- expt_info[[expt]]$meta_tps
  
  # Import clusters
  clusters <- read.table(sprintf("out/clusters_%s.txt", expt), header = T, sep = "\t", check.names = F)
  clusters <- sapply(1:nrow(clusters), function(i) strsplit(clusters$clones[i], ";")[[1]], simplify = F)
  cluster_sizes <- sapply(clusters, length)
  
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
  meta_counts <- read.table(sprintf("data/meta/meta_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  meta_depths <- read.table(sprintf("data/meta/meta_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
  meta_n_vars <- nrow(meta_counts)
  
  # Get diag vars in metagenomic data
  diag_vars_overlap <- geno_meta_overlap(diag_vars, geno_counts, meta_counts)
  geno_diag_ids <- diag_vars_overlap$geno_ids
  meta_diag_ids <- diag_vars_overlap$meta_ids
  
  # Get only clusters that have any variables at all
  included_clusters <- which(sapply(geno_diag_ids, function(x) !is.na(x[1])))
  
  # Subset metagenomic data
  subset_meta <- subset_meta_data(included_clusters, meta_counts, meta_depths, meta_diag_ids)
  
  # Count the number of timepoints with zero depth
  n_zeros <- sapply(subset_meta$meta_depths, function(x) sum(apply(x, 2, function(y) all(y == 0))))
  
  # Compute median depth
  median_depths <- sapply(subset_meta$meta_depths, median)
  
  # Count variants and produce output
  out <- data.frame(cluster = included_clusters,
                    cluster_size = cluster_sizes[included_clusters],
                    n_overlap_diag_vars = sapply(geno_diag_ids[included_clusters], length),
                    n_zero_depth_tps = n_zeros,
                    median_depth = median_depths)
  
  # Apply lineage filter
  out$keep <- out$n_zero_depth_tps <= max_zero_depths &
    out$n_overlap_diag_vars >= min_n_diag_vars & 
    out$median_depth >= min_n_diag_vars
  
  # Write inferred sequences
  write.table(out, sprintf("out/cluster_fit_info_%s.txt", expt), sep = "\t",
              row.names = F, quote = F)
}
