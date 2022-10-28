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
# Import full tree
full_tree <- read.tree("out/tree_all_clones_rerooted.tree")

# Import cluster colors
full_col_table <- read.table("out/col_table_full.txt", header = T)

# Plot circular tree
col_ploidies <- RColorBrewer::brewer.pal(3, "Set2")[1:2]
{
  pdf("out/tree_ploidy.pdf", width = 5.5, height = 5.5)
  ops <- par(mar = c(2, 2, 2, 2), xpd = T)
  temp_tree <- full_tree
  temp_tree$tip.label <- unlist(pretty_labels[temp_tree$tip.label])
  plot(temp_tree, type = "fan", cex = 0.5, label.offset = 0.03)
  add.scale.bar(-0.27, -0.26)
  legend(0.2, -0.19, legend = c("2N", "3N"), title = "Ploidy", col = 1,
         pch = 16, pt.cex = 1.75, bg = "transparent", box.col = "transparent")
  legend(0.2, -0.19, legend = c("2N", "3N"), title = "Ploidy", col = col_ploidies,
         pch = 16, pt.cex = 1.5, bg = "transparent", box.col = "transparent")
  tiplabels(bg = col_ploidies[(unlist(ploidies[full_tree$tip.label]) == 3) + 1],
            cex = 1.5, pch = 21, offset = 0.015)
  par(ops)
  dev.off()
}
