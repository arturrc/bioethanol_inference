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
  
  # Read count and depth data
  counts_dat <- read.table(sprintf("data/geno/geno_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  depths_dat <- read.table(sprintf("data/geno/geno_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
  n_vars <- nrow(counts_dat)
  
  # Import clusters
  clusters <- read.table(sprintf("out/clusters_%s.txt", expt), header = T, sep = "\t", check.names = F)
  clusters <- sapply(1:nrow(clusters), function(i) strsplit(clusters$clones[i], ";")[[1]], simplify = F)
  
  # Find diagnostic vars for all clusters
  diagnostic_vars <- sapply(clusters, function(cluster_clones){
    print(cluster_clones)
    
    # Check for ploidy consistency
    if(any(diff(unlist(ploidies[cluster_clones])) != 0)){
      return(NA)
    }
    
    # Partition data into cluster and other
    other_clones <- all_clones[!(all_clones %in% cluster_clones)]
    
    cluster_counts <- as.matrix(counts_dat[, cluster_clones])
    cluster_depths <- as.matrix(depths_dat[, cluster_clones])
    
    other_counts <- as.matrix(counts_dat[, other_clones])
    other_depths <- as.matrix(depths_dat[, other_clones])
    
    # Find diagnostic vars
    out <- find_diagnostic_vars(cluster_counts, cluster_depths, other_counts, other_depths,
                                ploidies[[cluster_clones[1]]], verbose = T, debug = F)
    
    if(length(out$diag_vars) == 0) out <- NA
    
    return(out)
  }, simplify = F)
  
  # Count number of diag vars
  n_diag_vars <- sapply(diagnostic_vars, function(x){
    if(is.na(x[1])) return(NA)
    return(length(x[[1]]))
  })
  
  # Create list of genotypes for all clusters
  genotypes <- sapply(seq_along(clusters), function(i){
    cluster_clones <- clusters[[i]][1]
    ploidy <- ploidies[[cluster_clones[1]]]
    out <- seq(0, 1, 1/ploidy)
    return(out)
  }, simplify = F)
  
  # Estimate genotype posterior probabilities with EM
  geno_log_posts <- sapply(seq_along(clusters), function(i){
    # Skip if no overlapping diag vars for cluster
    if(is.na(diagnostic_vars[[i]][i])) return(NA)
    
    # Subset geno and meta data
    cluster_clones <- clusters[[i]]
    geno_ids <- diagnostic_vars[[i]]$diag_vars
    cluster_geno_counts <- unname(apply(as.matrix(counts_dat[geno_ids, cluster_clones]), 1, sum))
    cluster_geno_depths <- unname(apply(as.matrix(depths_dat[geno_ids, cluster_clones]), 1, sum))
    
    # Flip genotype counts
    for(j in which(diagnostic_vars[[i]]$flip)){
      cluster_geno_counts[j] <- cluster_geno_depths[j] - cluster_geno_counts[j]
    }
    
    # Calculate geno posteriors for each mutation
    em1 <- em_geno(cluster_geno_counts, cluster_geno_depths, genotypes[[i]])
    out <- em1$log_posteriors
    return(out)
  }, simplify = F)

  # Remove diag vars that have highest genotype prob = 0
  diagnostic_vars <- sapply(seq_along(clusters), function(i){
    if(is.na(diagnostic_vars[i])) return(NA)
    
    keep_ids <- apply(geno_log_posts[[i]], 1, which.max) != 1
    if(sum(keep_ids) == 0) return(NA)
    
    out <- diagnostic_vars[[i]]
    out$diag_vars <- out$diag_vars[keep_ids]
    out$flip <- out$flip[keep_ids]
    out$heterogeneous_geno <- out$heterogeneous_geno[keep_ids]
    
    return(out)
  }, simplify = F)
  
  geno_log_posts <- sapply(seq_along(clusters), function(i){
    if(is.na(diagnostic_vars[i])) return(NA)
    
    keep_ids <- apply(geno_log_posts[[i]], 1, which.max) != 1
    if(sum(keep_ids) == 0) return(NA)
    
    out <- geno_log_posts[[i]][keep_ids,]
    if(!is.matrix(out)) out <- matrix(out, nrow = 1)

    return(out)
  }, simplify = F)
  
  # Prepare and write output
  included_clusters <- which(!is.na(diagnostic_vars))
  output <- sapply(included_clusters, function(cluster){
    out <- data.frame(cluster = cluster,
                      diag_var = diagnostic_vars[[cluster]]$diag_vars,
                      chrom = counts_dat$chrom[diagnostic_vars[[cluster]]$diag_vars],
                      pos = counts_dat$pos[diagnostic_vars[[cluster]]$diag_vars],
                      alt_allele = counts_dat$alt[diagnostic_vars[[cluster]]$diag_vars],
                      flip = diagnostic_vars[[cluster]]$flip,
                      geno_log_posts = apply(geno_log_posts[[cluster]], 1, function(x) paste0(x, collapse = ";")))
  }, simplify = F)
  output <- do.call(rbind, output)
  write.table(output, sprintf("out/diagnostic_vars_%s.txt", expt), sep = "\t",
              row.names = F, quote = F)
  
  # Prepare and write report
  report <- sapply(included_clusters, function(cluster){
    out <- data.frame(cluster = cluster,
                      n_diag_vars = length(diagnostic_vars[[cluster]]$diag_vars),
                      n_flip = sum(diagnostic_vars[[cluster]]$flip),
                      n_rejected_het = length(diagnostic_vars[[cluster]]$heterogeneous_geno),
                      n_dubious_geno = sum(apply(geno_log_posts[[cluster]], 1, function(x) sum(x > -2.3) > 1)))
  }, simplify = F)
  report <- do.call(rbind, report)
  write.table(report, sprintf("out/report_diagnostic_vars_%s.txt", expt), sep = "\t",
              row.names = F, quote = F)
}
