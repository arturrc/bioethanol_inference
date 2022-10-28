rm(list = ls())
options(stringsAsFactors = F)

# PREAMBLE =============================================
library(reshape2)
library(plyr)
library(ggplot2)
library(matrixStats)

source("vars.R")
source("functions.R")


rollapply <- function(y, x = NULL, window, FUN, simplify = T, na.rm = T){
  if(is.null(x)) x <- 1:(length(y) - window + 1)
  min_x <- min(x, na.rm = na.rm)
  max_x <- max(x, na.rm = na.rm)
  breaks <- seq(min_x, max_x, window)
  if(!(length(breaks) > 1)) stop("ERROR: window is too large")
  out <- sapply(2:length(breaks), function(i){
    if(i == 2){
      start_id <- 1
      end_id <- sum(x <= breaks[i])
    } else {
      start_id <- sum(x <= breaks[i - 1]) + 1
      end_id <- sum(x <= breaks[i])
    }
    return(FUN(y[start_id:end_id]))
  }, simplify = simplify)
  return(list(breaks = breaks, midpoints = (breaks[-length(breaks)] + breaks[-1])/2, z = out))
}

# CODE =================================================
# Create output directory
dir.create("out/ploidy", showWarnings = F)

# Define plotting variables
chr_lengths <- c(198713, 796422, 292921, 1523194, 565755, 270091, 1075965, 524445, 431614, 727129, 645009, 1059202, 917820, 778342, 1068661, 931290)
spacing <- 2e5
chr_lengths_cum_starts <- cumsum(c(0, chr_lengths[-16] + spacing))
chr_lengths_cum_ends <- chr_lengths_cum_starts + chr_lengths
chr_midpoins <- c(-spacing/2, chr_lengths_cum_ends + spacing/2)

for(expt in expts){
  # Import genotyping data
  counts_dat <- read.table(sprintf("data/geno/geno_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  depths_dat <- read.table(sprintf("data/geno/geno_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
  
  # Calculate allele frequencies
  af_dat <- counts_dat
  af_dat[,-(1:4)] <- af_dat[,-(1:4)]/depths_dat[,-(1:4)]
  
  # Define position along the genome for all variants
  cum_pos <- af_dat$pos + chr_lengths_cum_starts[af_dat$chrom]
  
  # Plot a single clone
  n_clones <- ncol(counts_dat) - 4
  for(clone_i in 1:n_clones){
    cat("\r", clone_i, "           ")
    clone_name <- names(af_dat)[4 + clone_i]
    pretty_clone_name <- pretty_labels[[clone_name]]
    file_clone_name <- pretty_clone_name
    file_clone_name <- gsub("*", "", file_clone_name, fixed = T)
    file_clone_name <- gsub(":", "-", file_clone_name, fixed = T)
    file_clone_name <- gsub("(", "-", file_clone_name, fixed = T)
    file_clone_name <- gsub(")", "", file_clone_name, fixed = T)
    png(sprintf("out/ploidy/%s.png", file_clone_name),
        width = 7, height = 4, units = "in", res = 200)
    ops <- par(mar = c(2, 2, 1, 1), adj = 0)
    layout(rbind(c(1, 3), c(2, 4)), widths = c(5, 1.5))
    
    plot(cum_pos, af_dat[,4 + clone_i], pch = 16, ylim = c(0, 1), cex = 0.2, col = rgb(0, 0, 0, 0.5),
         ylab = "AF", main = pretty_clone_name,
         xlab = "", axes = F)
    box()
    axis(2, at = c(0, 1/3, 1/2, 2/3, 1), labels = c("0", "", "", "", "1"))
    abline(v = chr_midpoins, col = "red")
    
    plot(cum_pos, depths_dat[,4 + clone_i], pch = 16, cex = 0.2, col = rgb(0, 0, 0, 0.5), ylim = c(0, median(depths_dat[,4 + clone_i])*3),
         ylab = "Coverage",
         xlab = "", axes = F)
    box()
    axis(2)
    abline(v = chr_midpoins, col = "red")
    
    cur_af <- af_dat[,4 + clone_i]
    cur_af <- cur_af[cur_af > 0 & cur_af < 1]
    hist(cur_af, breaks = seq(0, 1, 0.02), main = "", xlim = c(0.02, 0.98), col = 1, axes = F, xlab = "", ylab = "")
    axis(1, at = c(0, 1/3, 0.5, 2/3, 1), labels = c("0", "", "", "", "1"))
    axis(2)
    box()
    
    hist(depths_dat[,4 + clone_i], breaks = seq(0, 1000), xlim = c(0, 80), main = "", col = 1, axes = F, xlab = "", ylab = "")
    axis(1, at = c(0, 20, 40, 60, 80), labels = c("0", "", "40", "", "1"))
    axis(2)
    box()
    
    par(ops)
    dev.off()
  }
}
