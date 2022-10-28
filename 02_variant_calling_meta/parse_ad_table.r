rm(list = ls())
options(stringsAsFactors = F)

# PREAMBLE =============================================
args <- commandArgs(trailingOnly = T)
input_file <- args[1]
output_counts <- args[2]
output_depths <- args[3]

# CODE =================================================
# Import vcf
tabdat <- read.table(input_file, header = T, sep = "\t", check.names = F)
n_vars <- nrow(tabdat)

# Get sample names
sample_cols <- grep(".AD", names(tabdat))
sample_names <- substr(names(tabdat)[sample_cols], 1, nchar(names(tabdat)[sample_cols]) - 3)
n_samples <- length(sample_names)

# Create output tables
dummy <- as.data.frame(matrix(-1, nrow = nrow(tabdat), ncol = n_samples))
names(dummy) <- sample_names

counts_dat <- data.frame(chrom = tabdat$CHROM,
                         pos = tabdat$POS,
                         ref = tabdat$REF,
                         alt = tabdat$ALT,
                         dummy)
depths_dat <- counts_dat

# Filter out mutations with more than one variable site
keep_ids <- !grepl(",", counts_dat$alt)

# Parse ADs
for(i in which(keep_ids)){
  if(i%%100 == 0) cat(i, "\n")
  ads <- strsplit(unlist(tabdat[i, sample_cols]), ",")
  ads <- sapply(ads, as.numeric, simplify = F)
  counts_dat[i,-(1:4)] <- unname(sapply(ads, function(x) x[2]))
  depths_dat[i,-(1:4)] <- unname(sapply(ads, sum))
}

counts_dat <- counts_dat[keep_ids,]
depths_dat <- depths_dat[keep_ids,]

# Reorder timepoints
labels <- names(counts_dat)[-(1:4)]
out_counts <- data.frame(counts_dat[,1:4],
                         counts_dat[,order(as.numeric(substr(labels, 10, 100))) + 4])
out_depths <- data.frame(depths_dat[,1:4],
                         depths_dat[,order(as.numeric(substr(labels, 10, 100))) + 4])

# Rename chromosomes
chroms <- as.list(1:17)
names(chroms) <- c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                   "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                   "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                   "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4",
                   "NC_001224.1")
out_counts$chrom <- unname(unlist(chroms[out_counts$chrom]))
out_depths$chrom <- unname(unlist(chroms[out_depths$chrom]))

# Output counts and depths
write.table(out_counts, output_counts,
            sep = "\t", quote = F, row.names = F)
write.table(out_depths, output_depths,
            sep = "\t", quote = F, row.names = F)
