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
# Define function for calculating genotype posterior probability given 
# parameters and a uniform prior
calc_post_example <- function(ploidy, depth, e){
  genotypes <- (0:ploidy)/ploidy
  genotype_labels <- paste0(0:ploidy, "/3")
  priors <- rep(1/(ploidy + 1), ploidy + 1)
  
  lls <- matrix(NA, nrow = depth + 1, ncol = ploidy + 1)
  for(i in 1:nrow(lls)){
    a <- i - 1
    for(j in 1:ncol(lls)){
      g <- genotypes[j]
      lls[i, j] <- log_p_count_error(a, depth, g, e, n_sd = 1000)
    }
  }
  
  posts <- exp(t(apply(lls, 1, function(x) x - logSumExp(x))))
  
  plot_dat <- data.frame(counts = 0:depth, posts)
  names(plot_dat) <- c("counts", genotype_labels)
  plot_dat <- melt(plot_dat, id.vars = "counts", 
                   variable.name = "genotype", value.name = "post")
  plot_dat <- data.frame(ploidy = ploidy,
                         depth = depth,
                         e = e,
                         plot_dat)
  return(plot_dat)
}

# Calculate genotype posterior probabilities
make_cond <- function(ploidy, depth, e){
  data.frame(ploidy = ploidy, depth = depth, e = e)
}
conds <- make_cond(2, 20, 0)
conds <- rbind(conds, make_cond(2, 100, 0))
conds <- rbind(conds, make_cond(2, 20, 0.01))
conds <- rbind(conds, make_cond(2, 100, 0.01))
conds <- rbind(conds, make_cond(3, 20, 0))
conds <- rbind(conds, make_cond(3, 100, 0))
conds <- rbind(conds, make_cond(3, 20, 0.01))
conds <- rbind(conds, make_cond(3, 100, 0.01))

posts <- sapply(1:nrow(conds), function(i){
  cat(">", i, "\n")
  out <- calc_post_example(conds$ploidy[i], conds$depth[i], conds$e[i])
  return(out)
}, simplify = F)
posts <- do.call(rbind, posts)

posts$e <- factor(posts$e, levels = c(0.01, 0))
p1 <- ggplot(posts[posts$ploidy == 2,], aes(counts, post)) +
  facet_wrap(e~depth, scales = "free_x", 
             labeller = labeller(e = label_both, depth = label_both)) +
  geom_line(aes(color = genotype), size = 0.75) +
  scale_color_brewer(palette = "Set2") +
  scale_linetype_manual(breaks = c("0.01", "0"), values = c(2, 1)) +
  labs(x = "Alternate allele count", y = "Genotype posterior probability",
       color = "Genotype") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(hjust = 0, margin = margin(b = 2)))

ggsave("out/post_geno_prob_diploid.pdf", p1,
       width = 4, height = 4, units = "in")

p2 <- ggplot(posts[posts$ploidy == 3,], aes(counts, post)) +
  facet_wrap(e~depth, scales = "free_x", 
             labeller = labeller(e = label_both, depth = label_both)) +
  geom_line(aes(color = genotype), size = 0.75) +
  scale_color_brewer(palette = "Set2") +
  scale_linetype_manual(breaks = c("0.01", "0"), values = c(2, 1)) +
  labs(x = "Alternate allele count", y = "Genotype posterior probability",
       color = "Genotype") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(hjust = 0, margin = margin(b = 2)))

ggsave("out/post_geno_prob_triploid.pdf", p2,
       width = 4, height = 4, units = "in")
