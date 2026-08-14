# Stage B-5 plots: PCA/MDS/NJT — F-anchored cluster labels from
# data/processed/cluster_maf/cluster_labels.tsv (for the 501 ∩ subset).
# Lifted/simplified from scripts/gadi/ordination_plots.R.
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(ape); library(ggtree)
})
source(here::here("scripts/_setup.R"))

ORD <- file.path(DATA_PROC, "ordination")
labels <- read_tsv(file.path(DATA_PROC, "cluster_maf", "cluster_labels.tsv"),
                   show_col_types = FALSE) |>
  transmute(Sample = Sample, Cluster = Cluster)

# ---- PCA -------------------------------------------------------------------
pca <- read_table(file.path(ORD, "Pk.eigenvec"), col_names = FALSE, show_col_types = FALSE) |>
  select(-1) |>
  rename(Sample = X2) |>
  rename_with(~ paste0("PC", seq_along(.x)), starts_with("X"))
eigenval <- scan(file.path(ORD, "Pk.eigenval"), quiet = TRUE)
pve <- tibble(PC = seq_along(eigenval), pve = eigenval / sum(eigenval) * 100)

pca_lab <- pca |> left_join(labels, by = "Sample") |>
  mutate(Cluster = ifelse(is.na(Cluster), "Unlabelled (regen-only)", Cluster))

pca_plot <- ggplot(pca_lab, aes(PC1, PC2, colour = Cluster)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_colour_manual(values = c(cluster_cols, "Unlabelled (regen-only)" = "grey60")) +
  coord_equal() +
  xlab(sprintf("PC1 (%.1f%%)", pve$pve[1])) +
  ylab(sprintf("PC2 (%.1f%%)", pve$pve[2])) +
  ggtitle("PCA — F-anchored cluster labels (501 ∩ truth; 58 regen-only grey)") +
  theme_pk()
save_fig(pca_plot, "stage_b/13_pca_cluster", width = 8, height = 6)

pve_plot <- ggplot(pve[1:20, ], aes(PC, pve)) +
  geom_col(fill = "grey40") +
  ylab("% variance explained") + theme_pk()
save_fig(pve_plot, "stage_b/13_pve", width = 7, height = 4)

# ---- MDS -------------------------------------------------------------------
mds <- read_table(file.path(ORD, "Pk.mds"), col_names = TRUE, show_col_types = FALSE) |>
  transmute(Sample = IID, MDS1 = C1, MDS2 = C2) |>
  left_join(labels, by = "Sample") |>
  mutate(Cluster = ifelse(is.na(Cluster), "Unlabelled (regen-only)", Cluster))

mds_plot <- ggplot(mds, aes(MDS1, MDS2, colour = Cluster)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_colour_manual(values = c(cluster_cols, "Unlabelled (regen-only)" = "grey60")) +
  coord_equal() +
  ggtitle("MDS — F-anchored cluster labels") +
  theme_pk()
save_fig(mds_plot, "stage_b/13_mds_cluster", width = 8, height = 6)

# ---- NJT -------------------------------------------------------------------
dist_ids <- read_table(file.path(ORD, "Pk.dist.id"), col_names = FALSE, show_col_types = FALSE)$X1
dist_mat <- as.matrix(read.table(file.path(ORD, "Pk.dist")))
rownames(dist_mat) <- dist_ids
colnames(dist_mat) <- dist_ids
options(ignore.negative.edge = TRUE)
NJT <- ape::nj(as.dist(dist_mat))

NJT_meta <- tibble(Sample = dist_ids) |>
  left_join(labels, by = "Sample") |>
  mutate(Cluster = ifelse(is.na(Cluster), "Unlabelled (regen-only)", Cluster),
         Colour  = unname(c(cluster_cols, "Unlabelled (regen-only)" = "grey60")[Cluster]))

tree_plot <- ggtree(NJT, layout = "daylight", size = 0.3) %<+% NJT_meta +
  aes(colour = I(Colour))
save_fig(tree_plot, "stage_b/13_njt_cluster", width = 8, height = 8)

message(sprintf("[13b] PCA %d samples (of %d in .eigenvec) | MDS %d | NJT %d",
                nrow(pca_lab), nrow(pca), nrow(mds), length(dist_ids)))
message(sprintf("[13b] labelled samples on plots: %d (F-anchored 501 ∩)",
                sum(!is.na(labels$Cluster))))
