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
median_covs <- llply(expts, function(expt){
  depths_dat <- read.table(sprintf("data/geno/geno_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
  median_covs <- apply(depths_dat[,-c(1:4)], 2, median)
  out <- data.frame(expt = expt,
                    clone = names(depths_dat)[-(1:4)], 
                    median_cov = median_covs)
  out <- out[order(out$median_cov, decreasing = T),]
  out$clone <- factor(as.character(out$clone), levels = as.character(out$clone))
  row.names(out) <- NULL
  return(out)
})

range_median_covs <- range(sapply(median_covs, function(x) range(x$median_cov)))

for(i in 1:4){
  p <- ggplot(median_covs[[i]], aes(median_cov, clone)) +
    geom_point() +
    coord_cartesian(xlim = c(0, range_median_covs[2])) +
    theme_bw()
  ggsave(sprintf("out/median_coverage_rank_%s.pdf", expts[i]), p,
         width = 4, height = 5)
}



