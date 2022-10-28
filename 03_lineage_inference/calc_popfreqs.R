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
freq_ploidies_all <- list()
for(expt in expts){
  # Import clusters
  clusters <- read.table(sprintf("out/clusters_%s.txt", expt), header = T, sep = "\t", check.names = F)
  clusters <- sapply(1:nrow(clusters), function(i) strsplit(clusters$clones[i], ";")[[1]], simplify = F)
  cluster_sizes <- sapply(clusters, length)
  
  # Get cluster ploidies
  cluster_ploidies <- unname(unlist(ploidies[sapply(clusters, function(x) x[1])]))
  
  # Import inferred frequencies
  freqs <- read.table(sprintf("out/freqs_%s.txt", expt), header = T)
  included_clusters <- unique(freqs$cluster)
  
  # Get phylogenies from the data
  edges <- list_parents(clusters[included_clusters], included_clusters)
  
  # Find basal lineages, and calculate fraction population diploid
  # based on their inferred frequencies
  basal_lineages <- edges$child[is.na(edges$parent)]
  freqs_basal <- freqs[freqs$cluster %in% basal_lineages,]
  freqs_basal$ploidy <- cluster_ploidies[freqs_basal$cluster]
  freq_ploidies <- ddply(freqs_basal, .(ploidy, tp), function(x){
    data.frame(tp = x$tp[1],
               ploidy = x$ploidy[1],
               freq = sum(x$freq))
  })
  freq_ploidies <- ddply(freq_ploidies, "tp", function(x){
    x$freq_norm <- x$freq/sum(x$freq)
    mean_ploidy <- 1/sum(x$freq_norm/x$ploidy)
    x$popfreq <- mean_ploidy*x$freq_norm/x$ploidy
    return(x)
  })

  freq_ploidies_all <- c(freq_ploidies_all, list(data.frame(expt = expt, freq_ploidies)))
  
  # Compute population frequency of lineages
  popfreqs <- ddply(freqs, "tp", function(x){
    # browser()
    tp <- x$tp[1]
    freq_trip <- freq_ploidies$freq_norm[freq_ploidies$tp == tp & freq_ploidies$ploidy == 3]
    mean_ploidy <- 1/((1-freq_trip)/2 + freq_trip/3)
    x$freq <- mean_ploidy*x$freq/cluster_ploidies[x$cluster]
    return(x)
  })
  
  # Save pop freqs
  write.table(popfreqs, sprintf("out/popfreqs_%s.txt", expt), sep = "\t", quote = F, row.names = F)
}

freq_ploidies_all <- do.call(rbind, freq_ploidies_all)
site_names <- list("Site B - 2018", "Site B - 2019", "Site A - 2018", "Site A - 2019")
names(site_names) <- expts
freq_ploidies_all$site_name <- unlist(site_names[freq_ploidies_all$expt])
freq_ploidies_all$site_name <- factor(freq_ploidies_all$site_name, levels = site_names[c(3, 4, 1, 2)])
freq_ploidies_all$ploidy_pretty <- c("Diploid", "Triploid")[(freq_ploidies_all$ploidy == 3) + 1]

p <- ggplot(freq_ploidies_all, aes(tp)) +
  facet_wrap(site_name~., scales = "free_x") +
  geom_line(aes(x = tp, y = freq_norm, color = as.factor(ploidy_pretty), alpha = "metagenome")) +
  geom_line(aes(x = tp, y = popfreq, color = as.factor(ploidy_pretty), alpha = "population")) +
  labs(x = "Day", y = "Frequency", color = "Ploidy", alpha = "Frequency in") +
  scale_alpha_manual(breaks = c("metagenome", "population"),
                        values = c(0.25, 1)) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme(strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(hjust = 0))

ggsave("out/freq_ploidies.pdf", p,
       width = 4, height = 2.75, units = "in")
