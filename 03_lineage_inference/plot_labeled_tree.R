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
plot_inferred_only <- T
for(expt in expts){
  all_clones <- c(expt_info[[expt]]$starter_clones, expt_info[[expt]]$picked_clones)
  
  # Import rooted tree
  tree <- read.tree(sprintf("out/tree_%s.tree", expt))
  
  # Import color tables 
  col_table <- read.table(sprintf("out/col_table_%s.txt", expt), header = T)
  col_table$included <- T
  
  if(plot_inferred_only){
    # Import inferred frequencies
    freqs <- read.table(sprintf("out/freqs_%s.txt", expt), header = T)
    included_clusters <- unique(freqs$cluster)
    col_table$included <- col_table$cluster %in% included_clusters
  }
  
  # Plot labeled tree
  {
    temp_tree <- tree
    temp_tree$tip.label <- unlist(pretty_labels[temp_tree$tip.label])
    temp_tree$tip.label[substr(temp_tree$tip.label, 2, 2) == "1"] <- 
      substr(temp_tree$tip.label[substr(temp_tree$tip.label, 2, 2) == "1"], 5, 100)
    
    pdf(sprintf("out/tree_%s.pdf", expt), width = 3, height = 4)
    ops <- par(mar = c(0, 0, 0, 0), xpd = T)
    
    plot(temp_tree, cex = 0.65, label.offset = 0.015)
    
    add.scale.bar(cex = 0.5)
    tips <- 1:sum(col_table$is_tip)
    if(plot_inferred_only) tips <- tips[col_table$included[match(tips, col_table$node)]]
    tiplabels(tip = tips, 
              col = col_table$col[match(tips, col_table$node)],
              cex = 1.5, pch = 16, offset = 0.0075)
    nodes <- (sum(col_table$is_tip)+1):nrow(col_table)
    if(plot_inferred_only) nodes <- nodes[col_table$included[match(nodes, col_table$node)]]
    nodelabels(node = nodes,
               col = col_table$col[match(nodes, col_table$node)],
               cex = 1.5, pch = 16)
    
    par(ops)
    dev.off()
  }
}
