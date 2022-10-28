rm(list = ls())
options(stringsAsFactors = F)

# PREAMBLE =============================================
library(reshape2)
library(plyr)
library(matrixStats)
library(R.utils)

source("vars.R")
source("functions.R")

vars <- commandArgs(trailingOnly = T)
output <- vars[1]
inputs <- vars[-1]

# CODE =================================================
# Import data tables
dats <- sapply(inputs, function(x){
  read.table(x, header = T)
}, simplify = F)

# Merge tables
out <- do.call(rbind, dats)

# Save output
write.table(out, output, sep = "\t",
            quote = F, row.names = F)
