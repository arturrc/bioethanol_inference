# GENERAL USE ==========================================================================================
# Make sure object is a matrix
# Important for 1D data that sometimes comes as vectors.
# If making a vector into a matrix, set as either a row or column.
tomatrix <- function(x, vec_to_row = T){
  if(!is.matrix(x)){
    if(vec_to_row){
      x <- matrix(x, nrow = 1)
    } else {
      x <- matrix(x, ncol = 1)
    }
  }
  return(x)
}

# For a vector of p values, find the cutoff that lead to a given FOR
find_p_given_FOR <- function(ps, FOR = 0.05, p_cutoff = 0.5){
  n <- length(ps)
  n_false_est <- sum(ps > p_cutoff)/p_cutoff
  for(r in seq_along(ps)){
    curFOR <- 1 - (1 - ps[r])*n_false_est/(n - r + 1)
    if(curFOR < FOR) break
  }
  return(ps[r])
}


# SIMULATION =============================================================================================
# Simulate genotyping counts
sim_geno_counts <- function(n_mutations, genotype_depths, mutation_genotypes){
  out <- sapply(1:n_mutations, function(i){
    out <- rbinom(1, genotype_depths[i], mutation_genotypes[i])
    return(out)
  })
  return(out)
}

# Simulate metagenomics counts
sim_meta_counts <- function(mutation_genotypes, lineage_frequency, metagenomics_depth){
  out <- sapply(1:length(mutation_genotypes), function(i){
    out <- rbinom(1, metagenomics_depth, mutation_genotypes[i]*lineage_frequency)
    return(out)
  })
  return(out)
}

# LINEAGE INFERENCE ======================================================================================
# EM of genotyping data
log_p_count_dbinom <- function(geno_count, geno_depth, genotype){
  dbinom(geno_count, geno_depth, genotype, log = T)
}

log_p_count_error <- function(geno_count, geno_depth, genotype, e = 0.01, n_sd = 3){
  p <- 0.5
  n_min <- max(0, floor(geno_depth*e - n_sd*sqrt(geno_depth*e*(1-e))))
  n_max <- min(geno_depth, ceiling(geno_depth*e + n_sd*sqrt(geno_depth*e*(1-e))))
  out <- sapply(n_min:n_max, function(n){
    y_min <- max(0, floor(p*n - n_sd*sqrt(n*p*(1-p))))
    y_max <- min(n, geno_count, ceiling(p*n + n_sd*sqrt(n*p*(1-p))))
    out <- sapply(y_min:y_max, function(y){
      out <- dbinom(geno_count - y, geno_depth - n, genotype, log = T) +
        dbinom(y, n, p, log = T) +
        dbinom(n, geno_depth, e, log = T)
      return(out)
    })
    return(logSumExp(out))
  })
  return(logSumExp(out))
}

em_geno_iteration <- function(log_p_genotypes_estimate, geno_counts, geno_depths, genotypes, e = 0.01, n_sd = 3){
  U <- length(geno_counts)
  
  loglikelihoods <- matrix(NA, nrow = U, ncol = length(genotypes))
  for(m in 1:U){
    for(l in 1:length(genotypes)){
      if(e == 0){
        loglikelihoods[m, l] <- log_p_count_dbinom(geno_counts[m], geno_depths[m], genotypes[l])
      } else {
        loglikelihoods[m, l] <- log_p_count_error(geno_counts[m], geno_depths[m], genotypes[l], e = e, n_sd = n_sd)
      }
    }
  }
  
  logp_data <- apply(loglikelihoods, 1, function(x) logSumExp(x + log_p_genotypes_estimate))
  nll <- -sum(logp_data)
  
  logposteriors <- matrix(NA, nrow = U, ncol = length(genotypes))
  for(m in 1:U){
    for(l in 1:length(genotypes)){
      logposteriors[m, l] <- loglikelihoods[m,l] + log_p_genotypes_estimate[l] - logp_data[m]
    }
  }
  
  log_p_genotypes_estimate <- sapply(1:length(genotypes), function(i) logSumExp(logposteriors[,i]) - log(U))
  
  out <- list(log_p_genotypes_estimate = log_p_genotypes_estimate,
              log_posteriors = logposteriors,
              nll = nll)
  
  return(out)  
}

em_geno <- function(geno_counts, geno_depths, genotypes, e = 0.01, n_sd = 3, nll_threshold = 0.001, 
                    random_start = F, n_it = 20, debug = F, verbose = F){
  if(debug == T) browser()
  best_nll <- Inf
  out <- NULL
  
  if(!random_start) n_it <- 1
  for(it in 1:n_it){
    if(verbose) cat("> it = ", it, "\n", sep = "")
    if(random_start){
      p_genotypes_estimate <- c(0, runif(length(genotypes) - 1), 1) # Probability of each genotype
      p_genotypes_estimate <- diff(sort(p_genotypes_estimate))
    } else {
      p_genotypes_estimate <- rep(1/length(genotypes), length(genotypes))
    }
    log_p_genotypes_estimate <- log(p_genotypes_estimate)
    
    past_nll <- Inf
    while(T){
      cur_out <- em_geno_iteration(log_p_genotypes_estimate, geno_counts, geno_depths, genotypes, e = e, n_sd = n_sd)
      cur_nll <- cur_out$nll
      
      if(verbose) cat(c("nll = ", cur_nll, "\n"), sep = "")
      if(abs(past_nll - cur_nll) < nll_threshold) break
      
      log_p_genotypes_estimate <- cur_out$log_p_genotypes_estimate
      past_nll <- cur_nll
    }
    
    if(verbose) cat("=> p_gs = {", round(exp(cur_out$log_p_genotypes_estimate[1]), 3), 
                    ", ", round(exp(cur_out$p_genotypes_estimate[2]), 3), "}\n", sep = "")
    
    if(cur_nll < best_nll){
      if(verbose) cat("*Best model updated!\n", sep = "")
      best_nll <- cur_nll
      out <- cur_out
    }
  }
  return(out)
}

# Find most likely genotype for given count and depth data
most_likely_genotype <- function(count, depth, ploidy, e = 0.01, n_sd = 3, debug = F){
  if(debug) browser()
  genotypes <- (0:ploidy)/ploidy
  ps <- sapply(genotypes, function(g) log_p_count_error(count, depth, g, e = e, n_sd = n_sd))
  out <- genotypes[which.max(ps)]
  return(out)
}

# Test for genotype heterogeneity across a set of clones
# Given n clones, test whether at a given locus there is at least one of them
# with a genotype different from the others given counts and depth values at
# that locus.
test_genotype_heterogeneity <- function(counts, depths, ploidy, n_perms = 1000, exp_counts_thresh = 5){
  n_clones <- length(counts)
  
  joint_g <- most_likely_genotype(sum(counts), sum(depths), ploidy)
  
  exp_counts <- joint_g*depths
  x2 <- sum(((counts - exp_counts)^2)/exp_counts)
  
  sampling_vec <- sapply(1:n_clones, function(i) rep(i, depths[i]), simplify = F)
  sampling_vec <- do.call(c, sampling_vec)
  p_val <- NA
  if(all(exp_counts > exp_counts_thresh)){
    # do chi-square test
    p_val <- pchisq(x2, n_clones - 1, lower.tail = F)
  } else {
    # do permutation test
    simx2 <- sapply(1:n_perms, function(it){
      sampled_count_dests <- sample(sampling_vec, sum(counts), replace = F)
      sampled_counts <- rep(0, length(depths))
      for(i in seq_along(sampled_count_dests)){
        sampled_counts[sampled_count_dests[i]] <- sampled_counts[sampled_count_dests[i]] + 1
      }
      sampled_counts
      sampled_x2 <- sum(((sampled_counts - exp_counts)^2)/exp_counts)
      return(sampled_x2)
    })
    
    p_val <- sum(simx2 > x2)/n_perms
    
  }
  
  return(p_val)
}


# Find diagnostic variants for a given cluster given count data
find_diagnostic_vars <- function(cluster_counts, cluster_depths, other_counts = NULL, other_depths = NULL, cluster_ploidy,
                                 verbose = F, debug = F){
  if(debug) browser()
  # Validate object types
  cluster_counts <- tomatrix(cluster_counts, vec_to_row = F)
  cluster_depths <- tomatrix(cluster_depths, vec_to_row = F)

  # Allow compatibility when cluster includes all clones
  if(length(other_counts) == 0 | is.null(other_counts)){
    other_counts <- other_depths <- matrix(0, nrow = nrow(cluster_counts))
  }
  
  # Useful variables
  n_vars <- nrow(cluster_counts)
  cluster_size <- ncol(cluster_counts)
  
  # If cluster size is 1, then only check for counts
  if(cluster_size == 1){
    keep_muts <- sapply(1:n_vars, function(m){
      if((cluster_counts[m,] > 0) & all(other_counts[m,] == 0)){
        return(1)
      } else if((cluster_counts[m,] < cluster_depths[m,]) & all(other_counts[m,] == other_depths[m,])){
        return(2)
      } else {
        return(0)
      }
    })
    heterogeneous_geno <- NULL
    cluster_muts_ids <- which(keep_muts > 0)
    flip_bool <- keep_muts[cluster_muts_ids] == 2
    
  # Long method
  } else {
    cluster_muts <- rep(F, n_vars)
    p_vals <- rep(-1.1, n_vars)
    flip_bool <- rep(NA, n_vars)
    for(m in 1:n_vars){
      if(verbose){
        if(m %% 10000 == 1) cat(m, "\n")
      }
      
      # Is mutation unique to cluster and seen in all members of the cluster
      keep <- F
      if(all(cluster_counts[m,] > 0) & all(other_counts[m,] == 0)){
        flip_bool[m] <- F
        keep <- T
      } else if(all(cluster_counts[m,] < cluster_depths[m,]) & all(other_counts[m,] == other_depths[m,])){
        flip_bool[m] <- T
        keep <- T
      }
      
      if(keep){
        p_vals[m] <- test_genotype_heterogeneity(cluster_counts[m,], cluster_depths[m,], cluster_ploidy)
        cluster_muts[m] <- T
      }
    }
    
    p_vals_filt <- p_vals[cluster_muts]
    cluster_muts_filt <- which(cluster_muts)
    cluster_muts_filt <- cluster_muts_filt[order(p_vals_filt)]
    p_vals_filt <- sort(p_vals_filt)
    
    
    p_thresh <- find_p_given_FOR(p_vals_filt, FOR = 0.05)
    cluster_muts_ids <- sort(cluster_muts_filt[p_vals_filt >= p_thresh])
    heterogeneous_geno <- sort(cluster_muts_filt[p_vals_filt < p_thresh])
    flip_bool <- flip_bool[cluster_muts_ids]
  }
  
  out <- list(diag_vars = cluster_muts_ids,
              flip = flip_bool,
              heterogeneous_geno = heterogeneous_geno)
  
  return(out)
}

# Read diagnostic variants
read_diag_vars <- function(cluster_clones){
  filehandle <- paste0(sort(cluster_clones), collapse = "_")
  con <- file(sprintf("out/diagnostic_vars/%s-%s.txt", expt, filehandle), "r")
  diag_vars <- readLines(con)
  close(con)
  diag_vars <- as.integer(strsplit(diag_vars, "\t")[[1]])
  return(diag_vars)
}

# Likelihood model of lineage frequency
nll_lineage_freq <- function(lineage_freqs, meta_counts, meta_depths, geno_log_posts, genotypes){
  n_freqs <- length(lineage_freqs)
  lls <- rep(0, n_freqs)
  n_vars <- length(meta_counts)
  for(i in 1:n_freqs){
    for(m in 1:n_vars){
      curtest <- sapply(genotypes, function(gm) dbinom(meta_counts[m], meta_depths[m], gm*lineage_freqs[i], log = T))
      if(any(is.na(curtest))) stop("NAs produced by dbinom().")
      lls[i] <- lls[i] + 
        logSumExp(
          sapply(genotypes, function(gm) dbinom(meta_counts[m], meta_depths[m], gm*lineage_freqs[i], log = T)) + geno_log_posts[m,])
    }
  }
  
  return(-lls)
}

# Fit the lineage frequency to metagenomic data and genotype posteriors
fit_lineage_freqs <- function(meta_counts, meta_depths, geno_log_posts, genotypes){
  meta_counts <- tomatrix(meta_counts, vec_to_row = F)
  meta_depths <- tomatrix(meta_depths, vec_to_row = F)
  n_tps <- ncol(meta_counts)
  
  freqs <- sapply(1:n_tps, function(t){
    ids <- meta_depths[,t] > 0
    if(sum(ids) == 0) return(NA) # Return NA if no mutation with depth above 0
    
    cur_counts <- meta_counts[ids,t]
    cur_depths <- meta_depths[ids,t]
    
    out <- optim(0.5, function(freq) nll_lineage_freq(freq, cur_counts, cur_depths, geno_log_posts, genotypes),
          method = "Brent", lower = 0, upper = 1)$par
    return(out)
  })
  
  return(freqs)
}

# Infer frequencies of all lineages in the data independently
infer_freqs_ind <- function(meta_counts_table, meta_depths_table, meta_diag_ids, meta_diag_flip,
                            geno_log_posts, genotypes, 
                            timepoints, tp_ids = NULL,
                            max_zero_depths = Inf,
                            max_n_diag_vars = 5e2, random_seed = 42, verbose = F, debug = F){
  
  if(debug) browser()
  
  if(is.null(tp_ids)) tp_ids <- seq_along(timepoints)
  included_clusters <- which(sapply(meta_diag_ids, function(x) !is.na(list(x))))
  
  # Cycle through timepoints
  fits <- sapply(tp_ids, function(tp_id){
    if(verbose) cat("TIMEPOINT ", timepoints[tp_id], "\n", sep = "")
    
    # Cycle through clusters
    freqs <- sapply(included_clusters, function(i){
      # if(debug) browser()
      
      subset_meta <- subset_meta_data(i, meta_counts_table, meta_depths_table, meta_diag_ids, tp_id = tp_id)
      
      # Flip metagenome counts
      for(j in which(meta_diag_flip[[i]])){
        subset_meta$meta_counts[[1]][j] <- subset_meta$meta_depths[[1]][j] - subset_meta$meta_counts[[1]][j]
      }
      
      # Fit lineage frequencies
      out <- fit_lineage_freqs(subset_meta$meta_counts[[1]], subset_meta$meta_depths[[1]], geno_log_posts[[i]], genotypes[[i]])
      
      return(out)
    }, simplify = F)
    if(debug) browser()
    
    # Format output
    out <- data.frame(tp = timepoints[tp_id],
                      cluster = included_clusters,
                      freq = unlist(freqs))
    
    return(out)
  }, simplify = F)
  fits <- do.call(rbind, fits)
  return(fits)
}

# From a list of clusters, and their accompannying ids, output all direct child-parent relationships
list_parents <- function(clusters, cluster_ids){
  cluster_sizes <- sapply(clusters, length)
  parents <- sapply(seq_along(clusters), function(i){
    potential_ancestors <- which(cluster_sizes > cluster_sizes[i])
    if(length(potential_ancestors) == 0){
      return(NA)
    } else {
      ancestors <- sapply(potential_ancestors, function(j) all(clusters[[i]] %in% clusters[[j]]))
      ancestors <- potential_ancestors[ancestors]
      if(length(ancestors) == 0){
        return(NA)
      } else {
        parent <- ancestors[which.min(cluster_sizes[ancestors])]
        return(parent)
      }
    }
  })
  
  # Prepare NA parent output
  out <- data.frame(parent = NA, child = cluster_ids[is.na(parents)])
  
  # Prepare other output, if it exists
  if(nrow(out) < length(cluster_ids)){
    ids <- which(!is.na(parents))
    out2 <- data.frame(parent = parents[ids], 
                       child = cluster_ids[ids])
    out2 <- out2[order(sapply(out2$parent, function(i) cluster_sizes[i]), decreasing = T),]
    out2$parent <- cluster_ids[out2$parent] # change parent numbers to cluster ids
    out <- rbind(out, out2)
  }

  return(out)
}

# For a given dataframe of parent-child relationships, output all unique lines of ancestry
find_genealogies <- function(edges, na_parent = T){
  n_edges <- nrow(edges)
  covered_edges <- rep(F, n_edges)
  
  if(na_parent){
    root_rule <- function(x) is.na(x)
  } else {
    root_rule <- function(x) length(x) == 0
  }
  
  genealogies <- c()
  i <- 0
  go <- T
  while(go){
    i <- i + 1
    if(i > n_edges){
      go <- F
      next
    }
    if(covered_edges[i]) next
    
    covered_edges[i] <- T
    
    cur_genealogy <- edges[i,2]
    cur_parent <- edges[i,1]
    while(!root_rule(cur_parent)){
      cur_genealogy <- c(cur_genealogy, cur_parent)
      id <- which(edges[,2] == cur_parent)
      cur_parent <- edges[id, 1]
      covered_edges[id] <- T
    }
    
    genealogies <- c(genealogies, list(cur_genealogy))
  }
  
  # Filter out non-unique genealogies
  keep <- sapply(seq_along(genealogies), function(i){
    included <- sapply(seq_along(genealogies), function(j){
      all(genealogies[[i]] %in% genealogies[[j]])
    })
    if(sum(included) > 1){
      return(F)
    } else {
      return(T)
    }
  })
  
  return(genealogies[keep])
}

# Plot frequency data by genealogy
plot_genealogy_freqs <- function(freq_dat, cluster_info, genealogies, output_dir, create_dir = T, wipe_dir = F){
  if(!dir.exists(output_dir)){
    if(create_dir){
      dir.create(output_dir)
    } else {
      stop("Output directory doesn't exist!")
    }
  } else {
    if(wipe_dir){
      if(length(list.files(output_dir)) > 0) system2("rm", c("-r", sprintf("%s/*", output_dir)))
    }
  }
  
  plot_dat <- freq_dat
  cluster_sizes <- sapply(plot_dat$cluster, function(cluster_id) cluster_info$cluster_size[cluster_info$cluster == cluster_id])
  cluster_n_diag_vars <- sapply(plot_dat$cluster, function(cluster_id) cluster_info$n_diag_vars[cluster_info$cluster == cluster_id])
  plot_dat$label <- paste0(cluster_sizes, " / ", cluster_n_diag_vars)
  
  for(i in seq_along(genealogies)){
    genealogy <- genealogies[[i]]
    
    title <- sprintf("Clusters %s", paste0(genealogy, collapse = ", "))
    p <- ggplot(plot_dat[plot_dat$cluster %in% genealogy,], aes(tp, freq, group = cluster)) +
      geom_line(aes(color = label)) +
      labs(x = "Timepoint", y = "Inferred frequency", 
           color = "Cluster size/\n# of diag vars",
           title = title) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    file_path <- sprintf("%s/clusters_%s.pdf",
                         output_dir,
                         paste0(rev(genealogy), collapse = "_"))
    ggsave(file_path, p, width = 7, height = 3)
  }
}

# Heuristic for adjusting frequencies of lineages in relation
adjust_lineage_freqs <- function(freq_dat, freq_list, max_freqs, genealogies, debug = F){
  if(debug) browser()
  if(missing(freq_list))
    freq_list <- dlply(freq_dat$freqs, "cluster", function(x) return(x$freq[order(x$tp)]))
  
  # Adjust the frequencies of current sister clusters
  n_tps <- length(freq_list[[1]])
  
  if(length(max_freqs) == 1){
    max_freqs <- rep(max_freqs, n_tps)
  }
  
  sister_clusters <- unique(sapply(genealogies, function(x) tail(x, 1)))
  
  for(i in 1:n_tps){
    freqs <- sapply(sister_clusters, function(x) freq_list[[as.character(x)]][i])
    if(sum(freqs) > max_freqs[i]){
      freqs <- max_freqs[i]*freqs/sum(freqs)
      for(j in seq_along(sister_clusters)){
        freq_list[[as.character(sister_clusters[j])]][i] <- freqs[j]
      }
    }
  }
  
  # Now cycle through sister clusters and do the same at the lower level
  sister_clusters <- unique(sapply(genealogies, function(x) tail(x, 1)))
  freq_list <- sapply(sister_clusters, function(cluster){
    if(debug) browser()
    # Find genealogies starting on current cluster
    ids <- which(sapply(genealogies, function(x) cluster %in% x))
    
    # Define new genealogies, remove base cluster, remove empty lineages
    new_genealogies <- genealogies[ids]
    new_genealogies <- sapply(new_genealogies, function(x) head(x, -1), simplify = F)
    genealogy_lengths <- sapply(new_genealogies, length)
    new_genealogies <- new_genealogies[genealogy_lengths != 0]
    
    # If there are no lower levels, return current cluster frequency as final
    if(length(new_genealogies) == 0){
      out <- freq_list[as.character(cluster)] # returns a named list
      return(out)
    }
    
    # If there are lower levels, then adjust frequencies and return
    new_max_freqs <- freq_list[[as.character(cluster)]]
    new_clusters <- unique(unlist(new_genealogies))
    new_freq_list <- freq_list[as.character(new_clusters)]
    out <- adjust_lineage_freqs(freq_list = new_freq_list, max_freqs = new_max_freqs, genealogies = new_genealogies)
    out <- c(out, freq_list[as.character(cluster)])
    return(out)
  }, simplify = F)
  
  if(debug) browser()
  out_list <- c()
  for(x in freq_list){
    if(length(x) == 1) {
      out_list <- c(out_list, x[1])
    } else {
      temp <- sapply(1:length(x), function(i) x[[i]], simplify = F)
      names(temp) <- names(x)
      out_list <- c(out_list, temp)
    }
  }

  # Prepare output  
  if(missing(freq_dat)){
    return(out_list)
  } else {
    out <- freq_dat
    out$freqs <- ddply(out$freqs, "cluster", function(x){
      x <- x[order(x$tp),]
      x$freq <- out_list[[as.character(x$cluster[1])]]
      return(x)
    })
    return(out)
  }
}

# Find diagnostic variants in genotype data that also exist in the metagenome data
# Outputs indices of overlapping diagnostic variants in both genotype and metagenome data
geno_meta_overlap <- function(diag_vars, geno_counts, meta_counts, max_n_diag_vars = 5e2){
  diag_vars_only <- sapply(diag_vars, function(x){
    if(is.na(list(x))){
      return(NA)
    } else {
      return(x$diag_var)
    }
  })
  
  # Make geno and meta variant strings for comparison
  geno_var_ids <- sapply(1:nrow(geno_counts), function(i){
    paste0(c(geno_counts$chrom[i], geno_counts$pos[i], geno_counts$alt[i]), collapse = "_")
  })
  meta_var_ids <- sapply(1:nrow(meta_counts), function(i){
    paste0(c(meta_counts$chrom[i], meta_counts$pos[i], meta_counts$alt[i]), collapse = "_")
  })
  
  # Find out which clusters will be fit at all and their diagnostic snp indices
  diag_vars_overlap <- sapply(seq_along(diag_vars_only), function(i){
    # Setup empty output
    out <- list(geno_ids = NA, 
                meta_ids = NA)
    
    cluster_diag_vars <- diag_vars_only[[i]]
    if(is.na(cluster_diag_vars[1])) return(out) # skip if no diag vars
    
    # Find diagnostic vars present in metagenomic data
    overlap_var_ids <- mixedsort(intersect(meta_var_ids, geno_var_ids[cluster_diag_vars]))
    if(length(overlap_var_ids) == 0) return(out) # skip if no 
    
    # Limit number of diagnostic mutations
    if(length(overlap_var_ids) > max_n_diag_vars){
      overlap_var_ids <- mixedsort(sample(overlap_var_ids, max_n_diag_vars, replace = F))
    }
    
    # Get indices in both geno and meta data
    out$geno_ids <- sapply(overlap_var_ids, function(var_id) which(geno_var_ids == var_id), USE.NAMES = F)
    out$meta_ids <- sapply(overlap_var_ids, function(var_id) which(meta_var_ids == var_id), USE.NAMES = F)
    
    return(out)
  }, simplify = F)
  
  # Prepare output
  out_geno <- sapply(diag_vars_overlap, function(x) x$geno_ids, simplify = F)
  out_meta <- sapply(diag_vars_overlap, function(x) x$meta_ids, simplify = F)
  
  out <- list(geno_ids = out_geno,
              meta_ids = out_meta)
  
  return(out)
}

# Subset metagenome data
# Subset metagenome data to a set of clusters, and a set of diagnostic variants.
# By default, take all timepoints, but otherwise, you can also define the timepoint index that you
# would like to subset the data to.
subset_meta_data <- function(included_clusters, meta_counts_table, meta_depths_table, meta_diag_vars, tp_id = NULL){
  # Subset data
  subset_data <- sapply(seq_along(included_clusters), function(i){
    cluster_id <- included_clusters[i]
    
    # Subset geno and meta data
    if(is.null(tp_id)){
      cluster_meta_counts <- unname(as.matrix(meta_counts_table[meta_diag_vars[[cluster_id]], -(1:4)]))
      cluster_meta_depths <- unname(as.matrix(meta_depths_table[meta_diag_vars[[cluster_id]], -(1:4)]))
    } else {
      cluster_meta_counts <- unname(meta_counts_table[meta_diag_vars[[cluster_id]], 4 + tp_id])
      cluster_meta_depths <- unname(meta_depths_table[meta_diag_vars[[cluster_id]], 4 + tp_id])
    }
    
    # Prepare output
    out <- list(meta_counts = cluster_meta_counts,
                meta_depths = cluster_meta_depths)
    
    return(out)
  }, simplify = F)
  
  # Prepare output
  out_counts <- sapply(subset_data, function(x) x$meta_counts, simplify = F)
  out_depths <- sapply(subset_data, function(x) x$meta_depths, simplify = F)
  
  out <- list(meta_counts = out_counts,
              meta_depths = out_depths)
  
  return(out)  
}

# Draw random frequencies satisfying inequalities
draw_initial_freqs <- function(Umat, cvec, max_freq = 1, delta = 1e-3, verbose = F, debug = F, random_seed = NULL, max_it = 1e4){
  if(debug) browser()
  if(!is.null(random_seed)) set.seed(random_seed)
  n_clusters <- ncol(Umat)
  n_constraints <- ncol(Umat)
  rule <- function(guess){
    out <- all((Umat %*% guess - cvec) >= 0) & all(guess >= 0) & all(guess <= max_freq)
    return(out)
  }
  guess <- rep(0, n_clusters)
  
  # Make sure delta magnitude is not larger than max_freq
  if(delta > max_freq) delta <- max_freq
  
  # Move guess away from 0
  positive_it <- n_clusters*10
  i <- 0
  it <- 0
  while(i <= positive_it){
    it <- it + 1
    if(it > max_it) stop("Max number of iterations in drawing initial guess.")
    if(all(guess > 0)) i <- i + 1
    id <- floor(runif(1, 0, n_clusters)) + 1
    r <- runif(1, -delta, delta)
    new_guess <- guess
    new_guess[id] <- new_guess[id] + r
    if(rule(new_guess)){
      if(verbose) print(sprintf("OK: %s, %s", it, i))
      guess <- new_guess
    }
  }
  
  # Rescale guess randomly
  while(T){
    new_guess <- max_freq*runif(1)*guess/max(guess)
    if(rule(new_guess)){
      guess <- new_guess
      break
    }
  }
  
  return(guess)
}

# Joint inference function
# Optimizes frequencies jointly given parent-children inequality constraints in a single timepoint.
# meta_counts and meta_depths are lists with the metagenome data of a single timepoint over all included clusters.
infer_freqs_jointly <- function(meta_counts, meta_depths, geno_log_posts, genotypes, edges, cluster_ids, max_freq = 1,
                                n_its = 10, outer.eps = 1e-4, mu = 1e-5, grad_delta_f = 1e-5, min_freq = 1e-6,
                                constrOptim_timeout = 60, max_attempts = 10, nll_cutoff = 1,
                                random_seed = NULL, verbose = F, debug = F){
    
  if(debug) browser()
  
  # If max_freq < min_freq, output all freqs as max_freq
  # This step simplifies inference of very small lineages for which we don't care about anyways
  if(max_freq < min_freq){
    # Format and output data
    outfreq <- data.frame(cluster = cluster_ids,
                          freq = 0)
    
    outfit <- data.frame(outer.eps = outer.eps,
                         mu = mu,
                         n_its = NA,
                         n_attempts = NA,
                         comp_time = NA,
                         value = NA,
                         convergence = NA,
                         outer.iterations = NA,
                         barrier.value = NA)
    
    out <- list(freqs = outfreq, fit_stats = outfit)
    
    # Output NA freq for any cluster whose meta depths are all 0
    ids <- sapply(meta_depths, function(x) all(x == 0))
    out$freqs$freq[ids] <- NA
    
    return(out)
  }
  
  # Implement inequality constraints: Umat and cvec
  parents <- unique(edges$parent)
  
  n_constraints <- length(parents)
  n_clusters <- length(cluster_ids)
  
  Umat <- matrix(0, nrow = n_constraints, ncol = n_clusters)
  cvec <- rep(0, n_constraints)
  
  for(i in 1:n_constraints){
    parent <- parents[i]
    if(is.na(parent)){
      children <- edges$child[is.na(edges$parent)]
    } else {
      parent_id <- which(cluster_ids == parent)
      children <- edges$child[which(edges$parent == parent)]
    }
    children_ids <- which(cluster_ids %in% children)
    
    if(is.na(parent)){
      cvec[i] <- -max_freq
    } else {
      Umat[i,parent_id] <- 1
    }
    
    Umat[i,children_ids] <- -1
  }
  
  # Add positivity constraints to inequalities
  Umat2 <- matrix(0, nrow = n_clusters, ncol = n_clusters)
  diag(Umat2) <- 1
  cvec2 <- rep(0, n_clusters)
  
  Umat <- rbind(Umat, Umat2)
  cvec <- c(cvec, cvec2)
  
  # Define objective function
  f <- function(freqs){
    if(any(freqs < 0)) return(1e300) # Enforce positivity constraint; I get an error if I don't include this line
    
    nll <- 0
    for(i in seq_along(freqs)){
      cur_nll <- nll_lineage_freq(freqs[i], meta_counts[[i]], meta_depths[[i]], 
                                  geno_log_posts[[i]], genotypes[[i]])
      nll <- nll + cur_nll
    }
    
    return(nll)
  }
  
  # Define empirical gradient function
  g <- function(freqs, tp_id){
    while(any(freqs < grad_delta_f)) grad_delta_f <- grad_delta_f/2
    
    f_freqs <- f(freqs)
    out <- sapply(seq_along(freqs), function(i){
      cur_freqs <- freqs
      cur_freqs[i] <- cur_freqs[i] - grad_delta_f
      return((f_freqs - f(cur_freqs))/grad_delta_f)
    })
    
    return(out)
  }
  
  # Set RNG seed
  if(!is.null(random_seed)) set.seed(random_seed)
  
  # Iterate optimizer and find best solution
  best_obj <- Inf
  best_mlfit <- NULL
  best_i <- NA
  best_comptime <- NA
  for(it in 1:n_its){
    if(verbose) cat("> it ", it, "\n", sep = "")
    mlfit <- NULL
    i <- 0
    while(is.null(mlfit)){
      i <- i + 1
      if(verbose) cat("  attempt ", i, "\n", sep = "")
      
      if(i > max_attempts) stop("Maximum number of optimizing attempts reached.")
      
      withTimeout({
        initial_guess <- draw_initial_freqs(Umat, cvec, max_freq)
        try({
          time1 <- Sys.time()
          mlfit <- constrOptim(initial_guess, f, g, ui = Umat, ci = cvec, 
                               method = "BFGS", outer.eps = outer.eps, mu = mu, 
                               control = list(abstol = nll_cutoff))
          time2 <- Sys.time()
          comptime <- as.numeric(time2 - time1)
        })
      },timeout=constrOptim_timeout,
      onTimeout="warning")
    }
    if(mlfit$value < best_obj){
      if(verbose) cat("  --> new best NLL: ", mlfit$value, "\n", sep = "")
      best_obj <- mlfit$value
      best_mlfit <- mlfit
      best_i <- i
      best_comptime <- comptime
      if(best_obj < nll_cutoff) break
    }
  }
  
  # Format and output data
  outfreq <- data.frame(cluster = cluster_ids,
                        freq = best_mlfit$par)
  
  outfit <- data.frame(outer.eps = outer.eps,
                       mu = mu,
                       n_its = n_its,
                       n_attempts = best_i,
                       comp_time = best_comptime,
                       as.data.frame(best_mlfit[c("value", "convergence", "outer.iterations", "barrier.value")]))
  
  out <- list(freqs = outfreq, fit_stats = outfit)
  
  # Output NA freq for any cluster whose meta depths are all 0
  ids <- sapply(meta_depths, function(x) all(x == 0))
  out$freqs$freq[ids] <- NA
  
  # Bring all small frequencies to min_freq
  # This also deals with floating point issues that lead to small negative frequencies
  out$freqs$freq[out$freqs$freq < min_freq] <- 0
  
  return(out)
}

# Get children from a parent
get_children <- function(parent, edges){
  if(is.na(parent)){
    children <- sample(edges[,2][is.na(edges[,1])])
  } else {
    ids <- edges[,1] == parent
    ids[is.na(ids)] <- F
    children <- edges[,2][ids]
  }
  
  return(children)
}

# Get parent
get_parent <- function(child, edges){
  parent <- edges[,1][edges[,2] == child]
  return(parent)
}

# Get siblings
get_siblings <- function(child, edges){
  if(is.na(child)) return(integer())
  
  parent <- get_parent(child, edges)
  if(is.na(parent)){
    siblings <- edges[,2][is.na(edges[,1])]
  } else {
    ids <- edges[,1] == parent
    ids[is.na(ids)] <- F
    siblings <- edges[,2][ids]
  }
  siblings <- siblings[siblings != child]
  
  if(length(siblings) == 0) return(integer())
  
  return(siblings)
}

# Organize clusters in depth level
get_cluster_levels <- function(edges){
  genealogies <- find_genealogies(edges)
  genealogies_rev <- sapply(genealogies, rev, simplify = F)
  cluster_levels <- list()
  
  while(length(genealogies_rev) > 0){
    first_members <- sapply(genealogies_rev, head, n = 1)
    cur_level <- sort(unique(first_members))
    cluster_levels <- c(cluster_levels, list(cur_level))
    genealogies_rev <- sapply(genealogies_rev, function(x) x[-1], simplify = F)
    genealogy_lengths <- sapply(genealogies_rev, length)
    genealogies_rev <- genealogies_rev[genealogy_lengths > 0]
  }
  return(cluster_levels)
}

# Do joint inference recursively
# Start at highest lineage level and infer all parents and their children jointly.
infer_freqs_recursively <- function(meta_counts_table, meta_depths_table, meta_diag_ids, meta_flip,
                                    geno_log_posts, genotypes, edges, cluster_ids, 
                                    cluster_clones, timepoints, tp_ids = NULL,
                                    n_its = 10, outer.eps = 1e-4, mu = 1e-5, grad_delta_f = 1e-5, min_freq = 1e-6,
                                    constrOptim_timeout = 60, max_attempts = 10, nll_cutoff = 1,
                                    random_seed = NULL, verbose = F, debug = F){
  
  if(debug) browser()
  
  # Set RNG seed
  if(!is.null(random_seed)) set.seed(random_seed)
  
  # Cycle through timepoints
  if(is.null(tp_ids)) tp_ids <- seq_along(timepoints)
  fits <- sapply(tp_ids, function(tp_id){
    if(verbose) cat("TIMEPOINT ", timepoints[tp_id], "\n", sep = "")
    
    # Set empty objects with results
    inferred_freqs <- data.frame(tp = timepoints[tp_id],
                                 cluster = cluster_ids,
                                 freq = rep(NA, length(cluster_ids)))
    
    fit_stats <- list()
    
    # Choose parent fit order, randomizing among siblings
    genealogies <- find_genealogies(edges)
    genealogies_rev <- sapply(genealogies, rev, simplify = F)
    parents <- c(NA)
    
    while(T){
      genealogy_lengths <- sapply(genealogies_rev, length)
      genealogies_rev <- genealogies_rev[genealogy_lengths > 1]
      if(length(genealogies_rev) == 0) break
      first_members <- sapply(genealogies_rev, head, n = 1)
      cur_parents <- unique(first_members)
      if(length(cur_parents) == 1){
        parents <- c(parents, cur_parents)
      } else {
        parents <- c(parents, sample(cur_parents)) # Randomize order in which children will be fit
      }
      genealogies_rev <- sapply(genealogies_rev, function(x) x[-1], simplify = F)
    }
    
    # Do joint inference of all parents and their respective children
    for(parent in parents){
      if(!(parent %in% parents)) next # Skip if parent was ever removed from the list
      if(verbose) cat("# Parent ", parent, " ----------------------------------------------\n", sep = "")
      approved <- F
      while(!approved){
        # Define max frequency possible
        if(is.na(parent)){
          grandparent_room <- 1
        } else {
          grandparent <- edges$parent[edges$child == parent]
          if(is.na(grandparent)){
            grandparent_freq <- 1
            aunts <- edges$child[is.na(edges$parent)]
            aunts <- aunts[aunts != parent]
          } else {
            grandparent_freq <- inferred_freqs$freq[inferred_freqs$cluster == grandparent]
            ids <- edges$parent == grandparent
            ids[is.na(ids)] <- F
            aunts <- edges$child[ids]
            aunts <- aunts[aunts != parent]
          }
          if(length(aunts) == 0){
            grandparent_room <- grandparent_freq
          } else {
            aunt_freqs <- sapply(aunts, function(x) inferred_freqs$freq[inferred_freqs$cluster == x])
            grandparent_room <- grandparent_freq - sum(aunt_freqs, na.rm = T)
          }
        }
        
        # Define children
        if(is.na(parent)){
          children <- sample(edges$child[is.na(edges$parent)])
        } else {
          ids <- edges$parent == parent
          ids[is.na(ids)] <- F
          children <- edges$child[ids]
        }
        
        # Define current clusters
        if(is.na(parent)){
          cur_clusters <- children      
        } else {
          cur_clusters <- c(parent, children)
        }
        if(verbose) cat("doing clusters: ", cur_clusters, "\n")
        if(verbose) cat("max freq: ", grandparent_room, "\n")
        
        # Subset metagenomic data to included clusters' diag vars only
        subset_meta <- subset_meta_data(cur_clusters, meta_counts_table, meta_depths_table, meta_diag_ids, tp_id = tp_id)
        for(i in seq_along(cur_clusters)){
          cluster_id <- cur_clusters[i]
          for(mut_id in seq_along(subset_meta$meta_counts[[i]])){
            if(meta_flip[[cluster_id]][mut_id])
              subset_meta$meta_counts[[i]][mut_id] <- subset_meta$meta_depths[[i]][mut_id] - subset_meta$meta_counts[[i]][mut_id]
          }
        }
        
        # Subset edges to just current parent
        cur_edges <- list_parents(cluster_clones[cur_clusters], cur_clusters)
        
        # Do joint inference
        fits <- infer_freqs_jointly(subset_meta$meta_counts, subset_meta$meta_depths, 
                                    geno_log_posts[cur_clusters], genotypes[cur_clusters],
                                    cur_edges, cur_clusters, grandparent_room,
                                    n_its = n_its, constrOptim_timeout = constrOptim_timeout, 
                                    max_attempts = max_attempts, min_freq = min_freq,
                                    outer.eps = outer.eps, mu = mu, grad_delta_f = grad_delta_f,
                                    nll_cutoff = nll_cutoff,
                                    verbose = verbose, debug = debug, random_seed = NULL)
        
        # Verify if any lineages have NA freq, and if so, deal with it
        if(any(is.na(fits$freqs$freq))){
          # Is that lineage a parent? If so, eliminate from list of clusters and rerun inference including its children
          na_clusters <- fits$freqs$cluster[which(is.na(fits$freqs$freq))]
          parent_na_clusters <- na_clusters[na_clusters %in% edges$parent]
          if(length(parent_na_clusters) > 0){
            for(na_cluster in parent_na_clusters){
              cluster_ids <- cluster_ids[cluster_ids != na_cluster]
              parents <- parents[parents != na_cluster]
            }
            edges <- list_parents(cluster_clones[cluster_ids], cluster_ids)
            
            next
          } else {
            approved <- T
          }
        } else {
          approved <- T
        }
        
        # Update results
        for(i in seq_along(cur_clusters)){
          cluster_id <- cur_clusters[i]
          inferred_freqs$freq[inferred_freqs$cluster == cluster_id] <- fits$freqs$freq[i]
        }
        
        fit_stats <- c(fit_stats, list(data.frame(tp = timepoints[tp_id], parent = parent, fits$fit_stats)))
      }
    }
    
    # Prepare output
    fit_stats <- do.call(rbind, fit_stats)
    out <- list(freqs = inferred_freqs, fit_stats = fit_stats)
    
    return(out)
  }, simplify = F)
  
  # Prepare output
  freqs <- sapply(fits, function(x) x$freqs, simplify = F)
  freqs <- do.call(rbind, freqs)
  
  fit_stats <- sapply(fits, function(x) x$fit_stats, simplify = F)
  fit_stats <- do.call(rbind, fit_stats)
  
  out <- list(freqs = freqs, fit_stats = fit_stats)
  
  return(out)
}

# Plot all inferred frequencies in groups of parent and their children
plot_parent_children_freqs <- function(freqs, edges, cluster_info, output_dir, create_dir = T, wipe_dir = F){
  # Create or wipe directory
  if(!dir.exists(output_dir)){
    if(create_dir){
      dir.create(output_dir)
    } else {
      stop("Output directory doesn't exist!")
    }
  } else {
    if(wipe_dir){
      if(length(list.files(output_dir)) > 0) system2("rm", c("-r", sprintf("%s/*", output_dir)))
    }
  }
  
  # Get all parents in data
  parents <- unique(edges$parent)
  
  # Prepare plot data and label
  plot_dat <- freqs
  cluster_sizes <- sapply(plot_dat$cluster, function(cluster_id) cluster_info$cluster_size[cluster_info$cluster == cluster_id])
  cluster_n_diag_vars <- sapply(plot_dat$cluster, function(cluster_id) cluster_info$n_overlap_diag_vars[cluster_info$cluster == cluster_id])
  plot_dat$label <- paste0(plot_dat$cluster, " (", cluster_sizes, "; ", cluster_n_diag_vars, ")")
  plot_dat$keep <- sapply(plot_dat$cluster, function(cluster_id) cluster_info$keep[cluster_info$cluster == cluster_id])
  
  # Cycle through parents and plot
  for(i in seq_along(parents)){
    parent <- parents[i]
    children <- sort(get_children(parent, edges))
    
    cur_clusters <- c(parent, children)
    if(is.na(parent)){
      cols <- rep(RColorBrewer::brewer.pal(9, "Set1"), 10)[1:length(children)]
    } else {
      cols <- c(1, rep(RColorBrewer::brewer.pal(9, "Set1"), 10)[1:length(children)])
    }
    
    cur_labels <- sapply(cur_clusters, function(x) plot_dat$label[which(plot_dat$cluster == x)[1]])
    
    p <- ggplot(plot_dat[plot_dat$cluster %in% cur_clusters,], aes(tp, freq, group = cluster)) +
      geom_line(aes(color = label, size = label, lty = keep)) +
      guides(linetype = "none", size = "none") +
      labs(x = "Timepoint", y = "Inferred frequency", 
           color = "Cluster #\n(# of clones; # of muts)") +
      scale_color_manual(breaks = cur_labels, values = cols) +
      scale_size_manual(breaks = cur_labels, values = c(1.25, rep(0.75, length(cur_labels) - 1))) +
      scale_linetype_manual(breaks = c(T, F), values = c(1, 2)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    file_path <- sprintf("%s/clusters_%s.pdf",
                         output_dir,
                         paste0(cur_clusters, collapse = "_"))
    ggsave(file_path, p, width = 7, height = 3)
  }
}

# Validate frequencies
# Make sure that children inferred frequencies are smaller than parents'.
validate_frequencies <- function(freqs, edges, debug = F, tolerance = 1e-16){
  meta_tps <- sort(unique(freqs$tp))
  for(tp in meta_tps){
    x <- freqs[freqs$tp == tp,]
    parents <- unique(edges$parent)
    for(parent in parents){
      children <- get_children(parent, edges)
      children_freqs <- sapply(children, function(child) x$freq[x$cluster == child])
      if(is.na(parent)){
        parent_freq <- 1
      } else {
        parent_freq <- x$freq[x$cluster == parent]
      }
      
      if((sum(children_freqs) - parent_freq) > tolerance){
        if(debug) browser()
        return(F)
      }
    }
  }
  return(T)
}

# Find the tree node or tip ids given a list of clones (tips)
find_tree_ids <- function(rtree, clones_list){
  node_ids <- sapply(clones_list, function(x){
    ids <- which(rtree$tip.label %in% x)
    if(length(ids) == 1){
      out <- list(tree_id = ids,
                  node = F)
    } else {
      out <- list(tree_id = getMRCA(rtree, ids),
                  node = T)
    }
    return(out)
  }, simplify = F)
  node_ids <- unname(node_ids)
  return(node_ids)
}

# Get cluster plotting order that matches tree arrangement
# Requires lists of all clusters, and complete tree
get_plotting_order <- function(included_clusters, cluster_clones, rtree){
  # Get edges and genealogies
  edges <- list_parents(cluster_clones, seq_along(cluster_clones))
  genealogies <- find_genealogies(edges)
  genealogies_rev <- sapply(genealogies, rev, simplify = F)
  
  # Make tree id dictionary
  tree_ids <- sapply(find_tree_ids(rtree, clusters), function(x) x$tree_id)
  tree_ids <- as.list(tree_ids)
  names(tree_ids) <- as.character(seq_along(clusters))
  
  # Reorder genealogies based on tip tree id
  genealogies_rev <- genealogies_rev[order(sapply(genealogies_rev, function(x) tree_ids[[as.character(tail(x, 1))]]))]
  
  # Subset to included clusters
  genealogies_rev <- sapply(genealogies_rev, function(x) x[x %in% included_clusters], simplify = F)
  genealogies_rev <- genealogies_rev[sapply(genealogies_rev, length) != 0]
  
  # Get final plotting order
  out <- c()
  while(length(genealogies_rev) > 0){
    cur_clusters <- sapply(genealogies_rev, function(x) x[1])
    cur_clusters <- unique(cur_clusters)
    out <- c(out, list(cur_clusters))
    
    genealogies_rev <- sapply(genealogies_rev, function(x) x[-1], simplify = F)
    genealogies_rev <- genealogies_rev[sapply(genealogies_rev, length) != 0]
  }
  
  return(out)
}

# Get Muller plotting coordinates for a single timepoint
get_muller_posy <- function(freqs, edges, plotting_order, sibling_gap = T){
  ordered_lineages <- unlist(plotting_order)
  freqs$pos_y <- 0
  parents <- plotting_order[[1]]
  while(length(parents) > 0){
    # Make sure that parents follow correct order
    parents <- parents[order(sapply(parents, function(x) which(ordered_lineages == x)))]
    
    # Get current siblings to be positioned
    cur_lineages <- c(parents[1], get_siblings(parents[1], edges))

    # Make sure that cur_lineages follow correct order
    cur_lineages <- cur_lineages[order(sapply(cur_lineages, function(x) which(ordered_lineages == x)))]
        
    # Get current lineages' parent's position
    cur_parent <- get_parent(cur_lineages[1], edges)
    if(is.na(cur_parent)){
      start_pos_y <- 0
      parent_freq <- 1
    } else {
      id <- freqs$cluster == cur_parent
      start_pos_y <- freqs$pos_y[id]
      parent_freq <- freqs$freq[id]
    }
    
    # Calculate gap between cur_lineages
    cur_freqs <- sapply(cur_lineages, function(x){
      id <- freqs$cluster == x
      out <- freqs$freq[id]
      return(out)
    })
    
    if(sibling_gap){
      if(length(cur_lineages) > 1){
        freq_gap <- (parent_freq - sum(cur_freqs))/(length(cur_lineages) + 1)
      } else {
        freq_gap <- (parent_freq - sum(cur_freqs))/2
      }
      
    } else {
      freq_gap <- (parent_freq - sum(cur_freqs))/2
    }

    
    # Calculate lineages' positions
    cur_pos_ys <- c(start_pos_y + freq_gap)
    if(length(cur_lineages) > 1){
      for(i in 2:length(cur_lineages)){
        if(sibling_gap){
          new_freq <- tail(cur_pos_ys, 1) + cur_freqs[i - 1] + freq_gap
        } else {
          new_freq <- tail(cur_pos_ys, 1) + cur_freqs[i - 1]
        }
        cur_pos_ys <- c(cur_pos_ys, new_freq)
      }
    }
    
    # Add pos y to output dataframe
    for(i in seq_along(cur_lineages)){
      freqs$pos_y[freqs$cluster == cur_lineages[i]]<- cur_pos_ys[i]
    }
    
    # Update parents list
    parents <- parents[!(parents %in% cur_lineages)]
    for(i in seq_along(cur_lineages)){
      parents <- c(parents, get_children(cur_lineages[i], edges))
    }
  }
  
  return(freqs)
}

# Produce Muller plot
plot_muller <- function(freqs, edges, col_table, plotting_order, sibling_gap = T){
  # Get variables
  tps <- sort(unique(freqs$tp))
  ordered_clusters <- unlist(plotting_order)
  
  # Get pos_y for all timepoints
  posy <- ddply(freqs, "tp", function(x){
    get_muller_posy(x, edges, plotting_order, sibling_gap = sibling_gap)
  })
  
  # Make ggplot with polygons
  p <- ggplot()
  for(lineage in ordered_clusters){
    cur_freqs <- posy[posy$cluster == lineage,]
    cur_freqs <- cur_freqs[order(cur_freqs$tp),]
    xs <- c(tps, rev(tps))
    ys <- c(cur_freqs$pos_y, rev(cur_freqs$pos_y + cur_freqs$freq))
    poldat <- data.frame(cluster = lineage,
                         x = xs,
                         y = ys)
    p <- p + geom_polygon(data = poldat, aes(x = x, y = y, fill = as.factor(cluster)))
  }
  p <- p + 
    guides(fill = "none", color = "none")
  
  if(!is.null(col_table)){
    p <- p + scale_fill_manual(breaks = col_table$cluster, values = col_table$col)
  }
  
  return(p)
}

# Simplify and standardize clone names
simplify_tip_label <- function(x){
  if(grepl("PE-2", x) | grepl("SA-1", x) | grepl("FT-858", x) | grepl("IRA-D", x)){
    out <- paste0("*", x)
  } else if(grepl("_2019", x)){
    out <- substr(x, 1, regexpr("_", x) - 1)
  } else {
    out <- x
  }
  
  if(grepl("UIR", x) | grepl("UCP", x)){
    out <- substr(out, 4, 100)
  }
  
  out <- paste0(substr(out, 1, nchar(out) - 1), "(", substr(out, nchar(out), nchar(out)), ")")
  return(out)
}

