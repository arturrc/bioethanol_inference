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
# Import tree
tree <- read.tree("out/tree_all_clones_rerooted.tree")
n_tips <- length(tree$tip.label)
n_nodes <- tree$Nnode

# Assign clusters for the full tree
nodes <- c((n_tips + 1):(n_tips + n_nodes))
full_clusters <- sapply(nodes, function(node){
  extract.clade(tree, node)$tip.label  
}, simplify = F)
nodes <- c(nodes, 1:n_tips)
full_clusters <- c(full_clusters, as.list(tree$tip.label))

# Organize and name clusters
full_clusters <- sapply(full_clusters, function(x) mixedsort(x), simplify = F)
cluster_ids <- sapply(full_clusters, function(x) paste0(x, collapse = "_"))
ids <- mixedorder(cluster_ids)
full_clusters <- full_clusters[ids]
nodes <- nodes[ids]
ids <- order(sapply(full_clusters, length))
full_clusters <- full_clusters[ids]
nodes <- nodes[ids]
is_tip <- nodes <= n_tips

# Output clusters
output <- data.frame(cluster = seq_along(full_clusters),
                     clones = sapply(full_clusters, paste0, collapse = ";"),
                     node = nodes,
                     is_tip = is_tip)
write.table(output, "out/clusters_full.txt", sep = "\t",
            row.names = F, quote = F)

# Generate clusters for each experiment
for(expt in expts){
  # Import experiment tree
  tree <- read.tree(sprintf("out/tree_%s.tree", expt))
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  
  # Assign clusters for the full tree
  nodes <- (n_tips + 1):(n_tips + n_nodes)
  clusters <- sapply(nodes, function(node){
    extract.clade(tree, node)$tip.label  
  }, simplify = F)
  nodes <- c(nodes, 1:n_tips)
  clusters <- c(clusters, as.list(tree$tip.label))
  
  # Organize and name clusters
  clusters <- sapply(clusters, function(x) mixedsort(x), simplify = F)
  cluster_ids <- sapply(clusters, function(x) paste0(x, collapse = "_"))
  ids <- mixedorder(cluster_ids)
  clusters <- clusters[ids]
  nodes <- nodes[ids]
  ids <- order(sapply(clusters, length))
  clusters <- clusters[ids]
  nodes <- nodes[ids]
  
  # Find correspondence to full clusters
  match_full_clusters <- sapply(clusters, function(x){
    ids <- which(sapply(full_clusters, function(y){
      all(x %in% y)
    }))
    lengths <- sapply(full_clusters[ids], length)
    return(ids[which.min(lengths)])
  })
  
  # Output clusters
  output <- data.frame(cluster = seq_along(clusters),
                       clones = sapply(clusters, paste0, collapse = ";"),
                       node = nodes,
                       is_tip = nodes <= n_tips,
                       full_cluster = match_full_clusters)
  write.table(output, sprintf("out/clusters_%s.txt", expt), sep = "\t",
              row.names = F, quote = F)
}
