require('gtools')

# expt_info
expt_info_ <- read.table("data/info_expts.txt", header = T)
meta_info_ <- read.table("data/info_metagenomics.txt", header = T)
starter_strains_info_ <- read.table("data/info_starter_strains.txt", header = T)  
starter_clones_info_ <- read.table("data/info_starter_clones.txt", header = T)  
picked_clones_info_ <- read.table("data/info_picked_clones.txt", header = T)

# Set expt vector
expts <- expt_info_$expt

# ~> metadata info of each bioethanol production season
expt_info <- list()
for(expt_ in expts){
  meta_tps_ <- sort(meta_info_$timepoint[meta_info_$expt == expt_])
  out_starter_strains <- sort(starter_strains_info_$starter_strain[starter_strains_info_$expt == expt_])
  out_ <- list(site = expt_info_$site[expt_info_$expt == expt_],
              year = expt_info_$year[expt_info_$expt == expt_],
              meta_tps = meta_tps_,
              meta_samples = sapply(meta_tps_, function(tp_) meta_info_$label[meta_info_$expt == expt_ & meta_info_$timepoint == tp_]),
              starter_strains = out_starter_strains,
              starter_clones = mixedsort(starter_clones_info_$label[(starter_clones_info_$strain %in% out_starter_strains) & starter_clones_info_$keep]),
              starter_clones_rejected = mixedsort(starter_clones_info_$label[(starter_clones_info_$strain %in% out_starter_strains) & !starter_clones_info_$keep]),
              picked_clones = mixedsort(picked_clones_info_$label[(picked_clones_info_$expt == expt_) & picked_clones_info_$keep]),
              picked_clones_rejected = mixedsort(picked_clones_info_$label[(picked_clones_info_$expt == expt_) & !picked_clones_info_$keep]))
  
  expt_info[[expt_]] <- out_
}

# Get all included clones
all_included_clones <- sapply(expts, function(expt){
  c(expt_info[[expt]]$starter_clones, expt_info[[expt]]$picked_clones)
}, simplify = F)
all_included_clones <- unname(do.call(c, all_included_clones))
all_included_clones <- mixedsort(unique(all_included_clones))

# Get ploidies
ploidies <- as.list(c(picked_clones_info_$ploidy, starter_clones_info_$ploidy))
names(ploidies) <- c(picked_clones_info_$label, starter_clones_info_$label)

# Get pretty labels for clones
pretty_labels <- as.list(c(picked_clones_info_$pretty_label, starter_clones_info_$pretty_label))
names(pretty_labels) <- c(picked_clones_info_$label, starter_clones_info_$label)

# Remove temporary objects
rm(list = ls()[grep("_$", ls())])
