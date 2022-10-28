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
  all_clones <- c(expt_info[[expt]]$starter_clones, expt_info[[expt]]$picked_clones)
  tps <- expt_info[[expt]]$meta_tps
  
  # Read meta count and depth data
  meta_counts <- read.table(sprintf("data/meta/meta_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  meta_depths <- read.table(sprintf("data/meta/meta_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
  meta_n_vars <- nrow(meta_counts)
  
  # Calculate fraction matrix
  freq <- as.matrix(meta_counts[,-(1:4)]/meta_depths[,-(1:4)])
  
  # Produce plots
  png(sprintf("out/lines_%s.png", expt), width = 5, height = 3.5, units = "in", res = 900)
  par(mar = c(4.2, 4.2, 0.4, 0.4))
  ids <- sample(1:nrow(freq), 4000)
  plot(1, type = "n", xlim = c(0, max(tps))+c(max(tps)*0.03, -max(tps)*0.03), ylim = c(0.02, 0.98),
       xlab = "Day of sampling", ylab = "Frequency in metagenomic data")
  for(id in ids){
    lines(tps, freq[id,], lwd = 0.05, col = rgb(0, 0, 0, 0.7))
  }
  dev.off()
}
