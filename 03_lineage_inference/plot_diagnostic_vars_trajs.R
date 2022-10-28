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
  
  # Get cluster ploidies
  cluster_ploidies <- unname(unlist(ploidies[sapply(clusters, function(x) x[1])]))
  
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
  
  # Import inferred frequencies
  inf_freqs <- read.table(sprintf("out/freqs_%s.txt", expt), header = T)
  
  # Find diagnostic variant indices in both genotype and metagenome data
  diag_vars_overlap <- geno_meta_overlap(diag_vars, geno_counts, meta_counts)
  geno_diag_ids <- diag_vars_overlap$geno_ids
  meta_diag_ids <- diag_vars_overlap$meta_ids
  
  # Get set of included clusters
  included_clusters <- intersect(unique(inf_freqs$cluster),
                                 which(sapply(diag_vars, function(x) is.data.frame(x))))
  included_clusters <- intersect(included_clusters,
                                 which(sapply(geno_diag_ids, function(x) !is.na(x[1]))))
  
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
  
  # Get frequency trajectories
  freqs <- sapply(meta_diag_ids, function(x){
    if(is.na(x[1])) return(NA)
    out <- meta_counts[x, -(1:4)]/meta_depths[x, -(1:4)]
    return(out)
  }, simplify = F)
  
  # Get most likely genotype for each mutation
  diag_var_genos <- sapply(seq_along(clusters), function(i){
    if(is.na(geno_diag_ids[[i]][1])) return(NA)
    ids <- match(geno_diag_ids[[i]], diag_vars[[i]]$diag_var)
    out <- unname(apply(diag_vars[[i]][,-(1:6)], 1, which.max))
    return(out)
  }, simplify = F)
  
  # Flip trajectories
  for(i in included_clusters){
    for(j in which(diag_flip[[i]])){
      freqs[[i]][j,] <- 1 - freqs[[i]][j,]
    }
  }
  
  # Plot frequency trajectories
  genotype_cols <- c(rgb(1,0,0,0.7),
                     rgb(0,1,0,0.7),
                     rgb(0,0,1,0.7),
                     rgb(0.8,0.18,0.74,0.7))
  
  pdf(sprintf("out/diag_var_trajs_%s.pdf", expt), width = 7, height = 5)
  par(mfrow = c(1, 2))
  for(i in included_clusters){
    cols <- genotype_cols[diag_var_genos[[i]]]
    
    # plot in linear scale
    plot(1, type = "n", xlim = c(0, max(meta_tps)), ylim = c(0, 1),
         xlab = "Days", ylab = "Frequency in metagenomic data", 
         main = sprintf("%s: lineage %s (%sN)", expt, i, cluster_ploidies[i]))
    for(id in 1:nrow(freqs[[i]])){
      lines(meta_tps, freqs[[i]][id,], lwd = 0.5, col = cols[id])
    }
    ids <- which(inf_freqs$cluster == i)
    lines(inf_freqs$tp[ids], inf_freqs$freq[ids], lwd = 1.25, col = "grey20")
    
    # plot in log scale
    plot(1, type = "n", xlim = c(0, max(meta_tps)), ylim = c(1e-4, 1.01),
         xlab = "Days", ylab = "Frequency in metagenomic data", log = "y")
    for(id in 1:nrow(freqs[[i]])){
      lines(meta_tps, freqs[[i]][id,], lwd = 0.5, col = cols[id])
      points(meta_tps, freqs[[i]][id,], cex = 0.5, pch = 16, col = cols[id])
    }
    ids <- which(inf_freqs$cluster == i)
    lines(inf_freqs$tp[ids], inf_freqs$freq[ids], lwd = 1.25, col = "grey20")
  }
  dev.off()
}
