IBD.cutoffs <- c(0.125, 0.25, 0.5)

# include +-5% IBD
relaxed.IBD.cutoffs <- c(IBD.cutoffs, 1)
relaxed.IBD.cutoffs <- relaxed.IBD.cutoffs * 0.95




Pk.IBD.file <- "Mf.hmm_fract.txt"
Pk.IBD <- read.delim(Pk.IBD.file) 
  #mutate(sample1 = str_remove(sample1, "_DK.*")) %>% 
  #mutate(sample2 = str_remove(sample2, "_DK.*")) 
  

# Plot Fraction of sites that are IBD
IBD_meta_combined <- Pk.IBD %>%
    as.data.frame() %>%
    select(sample1, sample2, fract_sites_IBD) %>%
    left_join(
      metadata %>% mutate(sample1 = Sample) %>% select(sample1, Cluster, Proportion),
      by = "sample1"
    ) %>% 
    left_join(
      metadata %>% mutate(sample2 = Sample) %>% select(sample2, Cluster, Proportion),
      by = "sample2" 
    )
  
## Fraction IBD within and between Clusters
IBD_fract_plot <- IBD_meta_combined %>%
  unite("Clusters", c("Cluster.x", "Cluster.y"), sep = "-") %>% 
  ggplot(aes(x = Clusters, y = fract_sites_IBD, colour = Clusters)) +
  geom_boxplot() +
  theme(legend.position = "none", 
    legend.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  coord_flip() +
  scale_color_viridis_d() +
  ylab("Fraction of IBD sites")

ggsave("Mf_fract_IBD_clusters.png", dpi = 300, height = 7, width = 7, IBD_fract_plot)

# Base plots

Pk.label.cols <- c("Cluster")
Pk.legend.titles <- list("Cluster" = "Cluster")
Pk.prefixes <- metadata %>% 
  select(Cluster) %>%
  unique() %>% 
  as.list()


Pk.median.file <- "Mf_pairwise_IBD_median.txt"
Pk.median <- read.pairwise.matrix(Pk.median.file)

Pk.IBD.cutoffs <- c(fivenum(Pk.median), relaxed.IBD.cutoffs)

IBD.cutoffs <- c(0.002, 0.05, 0.01, .95)

plot.IBDs(
  Pk.prefixes,
  IBD.cutoffs,
  Pk.IBD,
  as.data.frame(metadata),
  Pk.label.cols,
  legend.titles = Pk.legend.titles
)

IBDs.plot.file <- "Mf_major_IBDs_cluster.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  as.data.frame(metadata),
  "Cluster",
  legend.title = "Cluster"
)



# Districts
metadata_dsitrict <- metadata %>% 
  mutate(sample1 = Sample) %>%
  mutate(sample1 = str_remove(sample1, "_DK.*")) %>% 
  mutate(sample1 = str_replace(sample1, "_", "-")) %>% 
  left_join(
    metadata_full %>% 
      select(Sample, district) %>% 
      rename(sample1 = Sample)
  ) %>%
  select(-sample1) %>% 
  rename(District = district)

Pk.label.cols <- c("District")
Pk.legend.titles <- list("District" = "District")
Pk.prefixes <- metadata_dsitrict %>% 
  select(District) %>%
  unique() %>% 
  as.list()

IBDs.plot.file <- "Mf_major_IBDs_dsitrict.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  as.data.frame(metadata_dsitrict),
  "District",
  legend.title = "District"
)




# Mn 
IBD.cutoffs <- c(0.125, 0.25, 0.5)

# include +-5% IBD
relaxed.IBD.cutoffs <- c(IBD.cutoffs, 1)
relaxed.IBD.cutoffs <- relaxed.IBD.cutoffs * 0.95

Pk.IBD.file <- "Mn.hmm_fract.txt"
Pk.IBD <- read.delim(Pk.IBD.file) 
  #mutate(sample1 = str_remove(sample1, "_DK.*")) %>% 
  #mutate(sample2 = str_remove(sample2, "_DK.*")) 
  

# Plot Fraction of sites that are IBD
IBD_meta_combined <- Pk.IBD %>%
    as.data.frame() %>%
    select(sample1, sample2, fract_sites_IBD) %>%
    left_join(
      metadata %>% mutate(sample1 = Sample) %>% select(sample1, Cluster, Proportion),
      by = "sample1"
    ) %>% 
    left_join(
      metadata %>% mutate(sample2 = Sample) %>% select(sample2, Cluster, Proportion),
      by = "sample2" 
    )
  
## Fraction IBD within and between Clusters
IBD_fract_plot <- IBD_meta_combined %>%
  unite("Clusters", c("Cluster.x", "Cluster.y"), sep = "-") %>% 
  ggplot(aes(x = Clusters, y = fract_sites_IBD, colour = Clusters)) +
  geom_boxplot() +
  theme(legend.position = "none", 
    legend.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  coord_flip() +
  scale_color_viridis_d() +
  ylab("Fraction of IBD sites")

ggsave("Mn_fract_IBD_clusters.png", dpi = 300, height = 7, width = 7, IBD_fract_plot)

# Base plots

Pk.label.cols <- c("Cluster")
Pk.legend.titles <- list("Cluster" = "Cluster")
Pk.prefixes <- metadata %>% 
  select(Cluster) %>%
  unique() %>% 
  as.list()


Pk.median.file <- "Mn_pairwise_IBD_median.txt"
Pk.median <- read.pairwise.matrix(Pk.median.file)

Pk.IBD.cutoffs <- c(fivenum(Pk.median), relaxed.IBD.cutoffs)

IBD.cutoffs <- c(0.002, 0.03, .95)

plot.IBDs(
  Pk.prefixes,
  IBD.cutoffs,
  Pk.IBD,
  as.data.frame(metadata),
  Pk.label.cols,
  legend.titles = Pk.legend.titles
)

IBDs.plot.file <- "Mn_major_IBDs_cluster.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  as.data.frame(metadata),
  "Cluster",
  legend.title = "Cluster"
)

# Districts
metadata_dsitrict <- metadata %>% 
  mutate(sample1 = Sample) %>%
  mutate(sample1 = str_remove(sample1, "_DK.*")) %>% 
  mutate(sample1 = str_replace(sample1, "_", "-")) %>% 
  left_join(
    metadata_full %>% 
      select(Sample, district) %>% 
      rename(sample1 = Sample)
  ) %>%
  select(-sample1) %>% 
  rename(District = district)

Pk.label.cols <- c("District")
Pk.legend.titles <- list("District" = "District")
Pk.prefixes <- metadata_dsitrict %>% 
  select(District) %>%
  unique() %>% 
  as.list()

IBDs.plot.file <- "Mn_major_IBDs_dsitrict.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  as.data.frame(metadata_dsitrict),
  "District",
  legend.title = "District"
)
