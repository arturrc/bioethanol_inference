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
  
  # Read geno count data
  geno_counts <- read.table(sprintf("data/geno/geno_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  geno_n_vars <- nrow(geno_counts)
  
  # Read meta count and depth data
  meta_counts <- read.table(sprintf("data/meta/meta_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  meta_depths <- read.table(sprintf("data/meta/meta_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
  meta_n_vars <- nrow(meta_counts)
  
  # Find diagnostic variant indices in both genotype and metagenome data
  diag_vars_overlap <- geno_meta_overlap(diag_vars, geno_counts, meta_counts)
  geno_diag_ids <- diag_vars_overlap$geno_ids
  meta_diag_ids <- diag_vars_overlap$meta_ids
  
  # Get boolean flip vectors for overlapping variants
  diag_flip <- sapply(seq_along(clusters), function(i){
    x <- geno_diag_ids[[i]]
    if(is.na(list(x))){
      return(NA)
    } else {
      out <- diag_vars[[i]]$flip[match(x, diag_vars[[i]]$diag_var)]
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
  
  # Create list of genotypes for all clusters
  genotypes <- sapply(seq_along(clusters), function(i){
    cluster_clones <- clusters[[i]][1]
    ploidy <- ploidies[[cluster_clones[1]]]
    out <- seq(0, 1, 1/ploidy)
    return(out)
  }, simplify = F)
  
  # Infer lineage frequencies
  freq_dat <- infer_freqs_ind(meta_counts, meta_depths, meta_diag_ids, diag_flip,
                              geno_log_posts, genotypes, 
                              meta_tps, tp_ids = NULL,
                              max_zero_depths = Inf,
                              max_n_diag_vars = 5e2, random_seed = 42, verbose = T, debug = F)
  
  # List all child-parent relationships
  edges <- list_parents(clusters[unique(freq_dat$cluster)], unique(freq_dat$cluster))
  
  # List all genealogies
  genealogies <- find_genealogies(edges)
  
  # Prepare cluster info
  cluster_info <- data.frame(cluster = unique(freq_dat$cluster),
                             cluster_size = sapply(unique(freq_dat$cluster), function(x) length(clusters[[x]])),
                             n_diag_vars = sapply(unique(freq_dat$cluster), function(x) length(meta_diag_ids[[x]])))
  
  # Plot inferred lineage frequencies
  plot_genealogy_freqs(freq_dat, cluster_info, genealogies, sprintf("out/genealogy_freqs_%s", expt),
                       create_dir = T, wipe_dir = T)
  
  # Write inferred sequences
  write.table(freq_dat, sprintf("out/inf_freqs_ind_%s.txt", expt), sep = "\t",
              row.names = F, quote = F)
}
