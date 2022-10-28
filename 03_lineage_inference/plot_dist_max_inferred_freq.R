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
freqs <- sapply(expts, function(expt){
  out <- read.table(sprintf("out/freqs_%s.txt", expt), header = T)
  out <- data.frame(expt = expt, out)
  return(out)
}, simplify = F)
freqs <- do.call(rbind, freqs)
row.names(freqs) <- NULL

length(unique(paste0(freqs$expt, "_", freqs$cluster)))

max_freqs <- ddply(freqs, .(expt, cluster), function(x){
  data.frame(expt = x$expt[1],
             cluster = x$cluster[1],
             max_freq = max(x$freq))
})

p <- ggplot(max_freqs, aes(max_freq)) +
  geom_histogram(bins = 20, fill = "grey20", color = "grey20") +
  labs(x = "Maximum inferred frequency", y = "No. of lineages") +
  scale_x_log10()

ggsave("out/hist_max_inferred_freq.pdf", p,
       width = 3, height = 2, units = "in")

