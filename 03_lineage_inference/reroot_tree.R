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
tree <- read.tree("data/tree1.ml.tree")

# Reroot tree at SA-1
rtree <- root(tree, node = 264)

# Subset to included clones only
rtree <- drop.tip(rtree, which(!(rtree$tip.label %in% all_included_clones)))

# Save tree
write.tree(rtree, "out/tree_all_clones_rerooted.tree")

# Subset experiment trees
for(expt in expts){
  all_clones <- c(expt_info[[expt]]$starter_clones, expt_info[[expt]]$picked_clones)
  
  # Remove any tips not in expt_info
  subtree <- drop.tip(rtree, which(!(rtree$tip.label %in% all_clones)))
  
  # Save tree
  write.tree(subtree, sprintf("out/tree_%s.tree", expt))
}

