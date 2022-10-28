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
# Get meta var ids
# Keeping only SNPs
meta_var_ids <- sapply(expts, function(expt){
  meta_counts <- read.table(sprintf("data/meta/meta_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  ids <- (meta_counts$ref %in% c("A", "T", "C", "G")) & (meta_counts$alt %in% c("A", "T", "C", "G"))
  meta_counts <- meta_counts[ids,]
  varids <- paste0(meta_counts$chrom, "_", meta_counts$pos, "_", meta_counts$ref, "_", meta_counts$alt)
  return(varids)
}, simplify = F)
names(meta_var_ids) <- expts

all_meta_var_ids <- unique(do.call(c, meta_var_ids))

# Get geno var ids
geno_var_ids <- sapply(expts, function(expt){
  counts_dat <- read.table(sprintf("data/geno/geno_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  varids <- paste0(counts_dat$chrom, "_", counts_dat$pos, "_", counts_dat$ref, "_", counts_dat$alt)
  return(varids)
}, simplify = F)
names(geno_var_ids) <- expts

all_geno_var_ids <- unique(do.call(c, geno_var_ids))

# Count SNPs and indels
is_snp_meta <- sapply(strsplit(all_meta_var_ids, "_"), function(x) nchar(x[3]) == 1 & nchar(x[4]) == 1)
is_snp_geno <- sapply(strsplit(all_geno_var_ids, "_"), function(x) nchar(x[3]) == 1 & nchar(x[4]) == 1)

sum(is_snp_meta)
sum(is_snp_geno)

# Count overlap between geno and meta
n_overlap <- sum(duplicated(c(all_geno_var_ids, all_meta_var_ids)))

# Count effective sequencing coverage
# Keeping only SNPs
meta_mean_depths <- sapply(expts, function(expt){
  meta_depths <- read.table(sprintf("data/meta/meta_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
  ids <- (meta_depths$ref %in% c("A", "T", "C", "G")) & (meta_depths$alt %in% c("A", "T", "C", "G"))
  meta_depths <- meta_depths[ids,]
  mean_depths <- apply(meta_depths[,-(1:4)], 2, mean, na.rm = T)
  return(mean_depths)
}, simplify = F)
sum(do.call(c, meta_mean_depths))

geno_mean_depths <- sapply(expts, function(expt){
  depths_dat <- read.table(sprintf("data/geno/geno_depths_%s.txt", expt), sep = "\t", header = T, check.names = F)
  mean_depths <- apply(depths_dat[,-(1:4)], 2, mean, na.rm = T)
  return(mean_depths)
}, simplify = F)
sum(do.call(c, geno_mean_depths))

# Count number of genes hit by these mutations
genome_annotation <- read.table("data/saccharomyces_cerevisiae_R64-3-1_20210421.gff", 
                                header = F, sep = "\t", quote = "",
                                skip = 21, nrows = 28386)
names(genome_annotation) <- c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")
genome_annotation$chr_n <- substr(genome_annotation$sequence, 4, 100)
genome_annotation$chr_n[genome_annotation$chr_n == "mt"] <- "XVII"
genome_annotation$chr_n <- as.numeric(as.roman(genome_annotation$chr_n))

yeast_genes <- genome_annotation[genome_annotation$feature == "gene",]
yeast_genes <- dlply(yeast_genes, "chr_n", function(x) x)

all_var_ids <- unique(c(all_geno_var_ids, all_meta_var_ids))

genes_hit <- sapply(yeast_genes, function(x) rep(F, nrow(x)), simplify = F)
for(x in strsplit(all_var_ids, "_")){
  chrom <- as.numeric(x[1])
  pos <- as.numeric(x[2])
  hits <- which(yeast_genes[[chrom]]$start <= pos & yeast_genes[[chrom]]$end >= pos)
  if(length(hits) > 0){
    genes_hit[[chrom]][hits] <- T
  }
}

sum(sapply(yeast_genes[-17], nrow))
sum(sapply(genes_hit[-17], sum))

# Count number of singletons
n_clones <- sapply(expts, function(expt){
  counts_dat <- read.table(sprintf("data/geno/geno_counts_%s.txt", expt), sep = "\t", header = T, check.names = F)
  varids <- paste0(counts_dat$chrom, "_", counts_dat$pos, "_", counts_dat$ref, "_", counts_dat$alt)
  n_clones <- apply(counts_dat[,expt_info[[expt]]$picked_clones], 1, function(x) sum(x > 0, na.rm = T))
  out <- data.frame(varids, n_clones)
  names(out) <- c("varid", sprintf("n_clone_%s", expt))
  return(out)
}, simplify = F)
n_clones_all <- merge(n_clones$UCP_2018, n_clones$UCP_2019, all = T)
n_clones_all <- merge(n_clones_all, n_clones$UIR_2018, all = T)
n_clones_all <- merge(n_clones_all, n_clones$UIR_2019, all = T)
n_clones_all[is.na(n_clones_all)] <- 0
n_clones_all$n_clones_total <- apply(n_clones_all[,-1], 1, sum)
n_clones_all <- n_clones_all[n_clones_all$n_clones_total != 0,]

sum(n_clones_all$n_clones_total == 1)
sum(n_clones_all$n_clones_total == max(n_clones_all$n_clones_total))/nrow(n_clones_all)

# Plot site frequency spectrum
max_n_clones <- max(n_clones_all$n_clones_total)
ids <- n_clones_all$n_clones_total > 1 & n_clones_all$n_clones_total < max_n_clones
hist(n_clones_all$n_clones_total[ids], breaks = max_n_clones*seq(0, 1, 1/30))

p <- ggplot(n_clones_all, aes(n_clones_total)) +
  stat_bin(breaks = seq(0, max_n_clones, 1), geom = "point", size = 0.5) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  labs(x = "No. of strains", y = "No. of variant sites") +
  scale_y_continuous(breaks = seq(0, 20000, 2000)) +
  scale_x_continuous(breaks = seq(0, 200, 20)) +
  coord_cartesian(expand = T, ylim = c(0, 16000), xlim = c(0, 135)) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("out/sfs_geno.pdf", p,
       width = 3, height = 2, units = "in")

# See overlap with 1011 var ids
varids1011 <- read.table("data/1011_varids.txt", header = T, sep = "\t")
varids1011$CHROM <- as.numeric(substr(varids1011$CHROM, 11, 12))
varids1011 <- with(varids1011, paste0(CHROM, "_", POS, "_", REF, "_", ALT))
ids <- match(all_var_ids, varids1011)
sum(!is.na(ids))

is_snp_1011 <- sapply(strsplit(varids1011, "_"), function(x) nchar(x[3]) == 1 & nchar(x[4]) == 1)
sum(is_snp_1011)
