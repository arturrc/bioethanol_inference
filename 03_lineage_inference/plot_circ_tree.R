rm(list = ls())
options(stringsAsFactors = F)

# PREAMBLE =============================================
library(reshape2)
library(plyr)
library(ggplot2)
library(matrixStats)
library(ape)

source("vars.R")
source("functions.R")

# CODE =================================================
node_numbers <- F
plot_inferred_only <- T

# Import full tree
full_tree <- read.tree("out/tree_all_clones_rerooted.tree")

# Import cluster colors
full_col_table <- read.table("out/col_table_full.txt", header = T)
full_col_table$included <- T

# Get inferred clusters
if(plot_inferred_only){
  included_clusters <- c()
  for(expt in expts){
    if(plot_inferred_only){
      # Import clusters
      clusters_dat <- read.table(sprintf("out/clusters_%s.txt", expt), header = T, sep = "\t", check.names = F)
      clusters <- sapply(1:nrow(clusters_dat), function(i) strsplit(clusters_dat$clones[i], ";")[[1]], simplify = F)
      cluster_sizes <- sapply(clusters, length)
      
      # Import inferred frequencies
      freqs <- read.table(sprintf("out/freqs_%s.txt", expt), header = T)
      
      # Get full cluster assignments of inferred clusters
      expt_clusters <- unique(freqs$cluster)
      included_clusters <- c(included_clusters, clusters_dat$full_cluster[match(expt_clusters, clusters_dat$cluster)])
    }
  }
  included_clusters <- unique(included_clusters)
  full_col_table$included <- full_col_table$cluster %in% included_clusters
}

# Plot circular tree
shape_ploidies <- c(5, 9)
{
  temp_tree <- full_tree
  temp_tree$tip.label <- unlist(pretty_labels[temp_tree$tip.label])
  
  pdf("out/tree_circ.pdf", width = 5.5, height = 6)
  ops <- par(mar = c(2, 2, 2, 2), xpd = T)
  
  plot(temp_tree, type = "fan", cex = 0.4, label.offset = 0.02, edge.color = "white")
  
  add.scale.bar(-0.27, -0.26)
  
  tips <- 1:sum(full_col_table$is_tip)
  if(plot_inferred_only) tips <- tips[full_col_table$included[match(tips, full_col_table$node)]]
  tiplabels(tip = tips, bg = full_col_table$col[match(tips, full_col_table$node)],
            cex = 1.5, pch = 21, offset = 0.01)
  
  nodes <- (sum(full_col_table$is_tip)+1):nrow(full_col_table)
  if(plot_inferred_only) nodes <- nodes[full_col_table$included[match(nodes, full_col_table$node)]]
  nodelabels(node = nodes, bg = full_col_table$col[match(nodes, full_col_table$node)],
            cex = 1.5, pch = 21)
  
  par(new = T)
  plot(temp_tree, type = "fan", cex = 0.4, label.offset = 0.02, edge.color = "white", edge.width = 2)
  par(new = T)
  plot(temp_tree, type = "fan", cex = 0.4, label.offset = 0.02, edge.color = "black", edge.width = 1)
  if(node_numbers){
    tiplabels(frame = "none", col = "white", cex = 0.3, adj = c(0.52, 0.5), offset = 0.015)
    tiplabels(frame = "none", col = "white", cex = 0.3, adj = c(0.48, 0.5), offset = 0.015)
    tiplabels(frame = "none", col = "white", cex = 0.3, adj = c(0.5, 0.53), offset = 0.015)
    tiplabels(frame = "none", col = "white", cex = 0.3, adj = c(0.5, 0.47), offset = 0.015)
    tiplabels(frame = "none", cex = 0.3, offset = 0.015)
    nodelabels(frame = "none", col = "white", cex = 0.3, adj = c(0.52, 0.5))
    nodelabels(frame = "none", col = "white", cex = 0.3, adj = c(0.48, 0.5))
    nodelabels(frame = "none", col = "white", cex = 0.3, adj = c(0.5, 0.53))
    nodelabels(frame = "none", col = "white", cex = 0.3, adj = c(0.5, 0.47))
    nodelabels(frame = "none", cex = 0.3)
  }

  # ploidy shapes
  ploidy_offset <- 0.023
  ploidy_offset <- 0.062
  ploidy_cex <- 0.85
  tiplabels(pch = 16 , col= "white", cex = 1, offset = ploidy_offset)
  tiplabels(pch = shape_ploidies[(unlist(ploidies[full_tree$tip.label]) == 3) + 1],
            cex = ploidy_cex, offset = ploidy_offset)
  legend(0.2, -0.19, legend = c("2N", "3N"), title = "Ploidy", col = 1,
         pch = shape_ploidies, pt.cex = ploidy_cex, bg = "transparent", box.col = "transparent")

  par(ops)
  dev.off()
}
