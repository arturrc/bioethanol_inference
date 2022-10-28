rm(list = ls())
options(stringsAsFactors = F)

# PREAMBLE =============================================
library(reshape2)
library(plyr)
library(ggplot2)
library(matrixStats)
library(ggmuller)
library(ape)

source("vars.R")
source("functions.R")
source("default_ggplot2_theme.R")

# CODE =================================================
node_rotations <- list("UCP_2018" = c(44, 45),
                       "UCP_2019" = c(),
                       "UIR_2018" = c(),
                       "UIR_2019" = c())
for(expt in expts){
  # Import clusters
  clusters <- read.table(sprintf("out/clusters_%s.txt", expt), header = T, sep = "\t", check.names = F)
  clusters <- sapply(1:nrow(clusters), function(i) strsplit(clusters$clones[i], ";")[[1]], simplify = F)
  cluster_sizes <- sapply(clusters, length)
  
  # Import inferred frequencies
  freqs <- read.table(sprintf("out/freqs_%s.txt", expt), header = T)
  included_clusters <- unique(freqs$cluster)
  
  # Import rooted tree
  rtree <- read.tree(sprintf("out/rooted_tree_%s.txt", expt))
  
  # Rotate tree nodes
  if(length(node_rotations[[expt]]) > 0){
    for(x in node_rotations[[expt]]){
      rtree <- rotate(rtree, x)
    }
    write.tree(rtree, file = "out/rtree.tmp")
    rtree <- read.tree("out/rtree.tmp")
  }
  
  # Import color tables
  col_table <- read.table(sprintf("out/col_table_%s.txt", expt), header = T)
  
  # Get phylogenies from the data
  edges <- list_parents(clusters[included_clusters], included_clusters)
  
  # Determine plotting order
  plotting_order <- get_plotting_order(included_clusters, clusters, rtree)
  
  # Validate frequencies
  if(!validate_frequencies(freqs, edges)) stop("Frequencies not valid!")
  
  # Produce Muller ggplot
  p <- plot_muller(freqs, edges, col_table, plotting_order, sibling_gap = F)
  
  p <- p +
    labs(x = "Day of sampling", y = "Frequency in the population") +
    scale_x_continuous(limits = c(1, max(freqs$tp)),
                       breaks = c(1, pretty(c(0, max(freqs$tp))))) +
    scale_y_continuous(limits = c(0, 1)) +
    coord_cartesian(expand = F, clip = "off") +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey90", colour = "black"))
  
  ggsave(sprintf("out/muller_%s.pdf", expt), p,
         width = 3.65, height = 2.3)
}
