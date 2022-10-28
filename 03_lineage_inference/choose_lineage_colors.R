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

# FUNCTIONS ===========================================
# Color conversion
hex_to_hsv <- function(hex_col) rgb2hsv(col2rgb(hex_col))
hsv_to_hex <- function(hsv_col) hsv(hsv_col[1], hsv_col[2], hsv_col[3],)

# Darken or lighten color
# factor > 1 -> lighten
alt_value_ <- function(hex_col, factor = 1){
  hsv_col <- hex_to_hsv(hex_col)
  if(factor > 1){
    hsv_col[3] <- hsv_col[3] + (1 - hsv_col[3])*(factor - 1)
  } else {
    hsv_col[3] <- hsv_col[3]*factor
  }
  return(hsv_to_hex(hsv_col))
}
alt_value <- Vectorize(alt_value_, USE.NAMES = F)

# Make color more or less saturated
# factor > 1 -> saturate
alt_sat_ <- function(hex_col, factor = 1){
  hsv_col <- hex_to_hsv(hex_col)
  if(factor > 1){
    hsv_col[2] <- hsv_col[2] + (1 - hsv_col[2])*(factor - 1)
  } else {
    hsv_col[2] <- hsv_col[2]*factor
  }
  return(hsv_to_hex(hsv_col))
}
alt_sat <- Vectorize(alt_sat_, USE.NAMES = F)


# Change hue
# Move hue by a certain number of degrees in the hue wheel (where 360 is a full circle)
alt_hue <- function(hex_col, degrees = 1){
  hsv_col <- hex_to_hsv(hex_col)
  hsv_col[1] <- hsv_col[1] + degrees/360
  hsv_col[1] <- hsv_col[1]%%1
  return(hsv_to_hex(hsv_col))
}

# Vary hue randomly
# Samples normally-distributed hues around a central hue
var_hue_norm <- function(hex_col, n = 1, sd = 10, seed = NULL){ # sd in degrees
  if(!is.null(seed)) set.seed(seed)
  sapply(1:n, function(i) alt_hue(hex_col, rnorm(1, 0, sd)))
}

# Get all tips from a node
tips_from_node <- function(tree, node) extract.clade(tree, node)$tip.label

# Get all subnodes and tips from given node
subnodes <- function(tree, node){
  tips <- tips_from_node(tree, node)
  nodes <- sapply(tips, function(x){
    nodepath(tree, node, match(x, tree$tip.label))[-1]
  }, simplify = F)
  nodes <- do.call(c, nodes)
  nodes <- sort(unique(nodes))
  return(nodes)
}

# CODE =================================================
# Import full tree
tree <- read.tree("out/tree_all_clones_rerooted.tree")
n_tips <- length(tree$tip.label)

# Import full clusters
full_clusters_dat <- read.table("out/clusters_full.txt", header = T, sep = "\t", check.names = F)
full_clusters <- sapply(1:nrow(full_clusters_dat), function(i) strsplit(full_clusters_dat$clones[i], ";")[[1]], simplify = F)

# Set empty color table
na_col <- "grey80"
full_col_table <- data.frame(cluster = seq_along(full_clusters),
                        node = full_clusters_dat$node,
                        is_tip = full_clusters_dat$is_tip,
                        col = na_col)

# Set basic color vectors
basic_cols <- c("brown3",
                RColorBrewer::brewer.pal(8, "Set1")[-1], 
                RColorBrewer::brewer.pal(6, "Dark2"), 
                RColorBrewer::brewer.pal(8, "Accent")[5:7],
                "chartreuse4")

plot(seq_along(basic_cols), rep(1, length(basic_cols)), pch = 16, cex = 4, col = basic_cols)

# # Set some basic color scheme
# full_col_table <- full_col_table[order(full_col_table$node),]
# full_col_table$col[!full_col_table$is_tip] <- rainbow(sum(!full_col_table$is_tip), s = 0.8, v = 0.9)
# full_col_table$col[full_col_table$is_tip] <- sapply(which(full_col_table$is_tip), function(i){
#   full_col_table$col[full_col_table$node == nodepath(tree, full_col_table$node[i], n_tips + 1)[2]]
# })


# Set colors
node_color <- function(node, full_col_table){
  full_col_table$col[full_col_table$node == node]
}

color_node <- function(node, col, full_col_table){
  full_col_table$col[match(node, full_col_table$node)] <- col
  return(full_col_table)
}

color_subtree <- function(node, base_col, sd, value_factor = 1, tree, full_color_table, seed){
  full_col_table$col[full_col_table$node == node] <- base_col
  nodes <- subnodes(tree, node)
  cols <- var_hue_norm(base_col, length(nodes), sd, seed)
  cols <- alt_value(cols, factor = value_factor)
  full_col_table$col[match(nodes, full_col_table$node)] <- cols
  return(full_col_table)
}

full_col_table <- color_subtree(152, alt_value(basic_cols[17], 1.1), 15, 0.9, tree, full_col_table, 2) # FT-858
full_col_table <- color_subtree(153, alt_value(node_color(152, full_col_table), 0.95), 15, 0.8, tree, full_col_table, 3)

full_col_table$col[full_col_table$node == 14] <- alt_value(basic_cols[1], 0.85) # A19-134-a # 1 but darker

full_col_table <- color_node(165, alt_sat(alt_value(basic_cols[1], 1.05), 1.05), full_col_table) # PE-2
full_col_table <- color_subtree(186, var_hue_norm(node_color(165, full_col_table), 1, 30, 101), 8, 0.9, tree, full_col_table, 25) # PE-2 subtree 1
full_col_table <- color_subtree(166, var_hue_norm(node_color(165, full_col_table), 1, 30, 102), 8, 0.9, tree, full_col_table, 25) # PE-2 subtree 2

full_col_table <- color_subtree(210, alt_value(basic_cols[4], 1.2), 15, 0.8, tree, full_col_table, 101) # IRA-D
full_col_table <- color_subtree(211, alt_value(node_color(210, full_col_table), 0.96), 15, 0.8, tree, full_col_table, 2134)
full_col_table <- color_subtree(212, alt_value(node_color(211, full_col_table), 0.96), 15, 0.8, tree, full_col_table, 454)
full_col_table <- color_subtree(213, alt_value(node_color(212, full_col_table), 0.96), 15, 0.8, tree, full_col_table, 13)
full_col_table <- color_subtree(214, alt_value(node_color(213, full_col_table), 0.96), 15, 0.8, tree, full_col_table, 2345)

full_col_table <- color_subtree(237, basic_cols[14], 8, 0.9, tree, full_col_table, 33) # Invaders 1

full_col_table <- color_node(244, alt_value(basic_cols[15], 1), full_col_table)
full_col_table <- color_node(249, alt_value(node_color(244, full_col_table), 0.97), full_col_table)
full_col_table <- color_node(259, alt_value(node_color(249, full_col_table), 0.97), full_col_table)

full_col_table <- color_subtree(245, basic_cols[4], 15, 0.9, tree, full_col_table, 31) # Invaders 2

full_col_table <- color_node(100, alt_hue(alt_value(basic_cols[4], 0.85), 45), full_col_table)

full_col_table <- color_subtree(250, alt_value(basic_cols[2], 0.9), 15, 0.9, tree, full_col_table, 31) # Invaders 3

full_col_table <- color_subtree(260, alt_sat(basic_cols[12], 0.85), 15, 0.9, tree, full_col_table, 32) # Invaders 4

full_col_table <- color_subtree(261, alt_value(basic_cols[8], 0.9), 15, 0.9, tree, full_col_table, 29)  # Invaders 5

full_col_table <- color_subtree(265, alt_value(basic_cols[7], 1.5), 15, 0.9, tree, full_col_table, 34) # Invaders 6

full_col_table <- color_node(123, alt_value(basic_cols[5], 0.85), full_col_table) # B18-87-a

full_col_table <- color_subtree(269, basic_cols[3], 15, 0.85, tree, full_col_table, 35) # Unknown

full_col_table <- color_subtree(284, alt_value(basic_cols[9], 1.1), 15, 0.9, tree, full_col_table, 36) # SA-1

# Save full color table
write.table(full_col_table, "out/col_table_full.txt", sep = "\t",
            quote = T, row.names = F)

# Save individual experiment color tables
for(expt in expts){
  # Import full clusters
  clusters_dat <- read.table(sprintf("out/clusters_%s.txt", expt), header = T, sep = "\t", check.names = F)
  clusters <- sapply(1:nrow(clusters_dat), function(i) strsplit(clusters_dat$clones[i], ";")[[1]], simplify = F)
  
  # Make col_table
  col_table <- data.frame(cluster = seq_along(clusters),
                          node = clusters_dat$node,
                          is_tip = clusters_dat$is_tip,
                          col = na_col)
  
  col_table$col <- sapply(clusters_dat$full_cluster, function(i) full_col_table$col[full_col_table$cluster == i])
  
  # Save expt color table
  write.table(col_table, sprintf("out/col_table_%s.txt", expt), sep = "\t",
              quote = T, row.names = F)
}

