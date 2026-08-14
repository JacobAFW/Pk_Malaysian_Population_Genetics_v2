# Stage B-6 supplement: plot connected components in IBD_Clonal_transmission_miss5.tsv
# (the truth's clonal-pair table) to show clonal clusters per macaque host.
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(igraph)
})
source(here::here("scripts/_setup.R"))

clonal <- read_tsv(file.path(DATA_RAW, "hmmIBD", "IBD_Clonal_transmission_miss5.tsv"),
                   show_col_types = FALSE)

cat("clonal-transmission pairs by tag:\n")
print(table(clonal$clonal_cluster))

for (tag in unique(clonal$clonal_cluster)) {
  sub <- clonal |> filter(clonal_cluster == tag)
  g <- graph_from_data_frame(sub |> transmute(from = sample1, to = sample2, weight = fract_sites_IBD),
                             directed = FALSE)
  comp <- components(g)
  n_multi <- sum(comp$csize > 1)
  V(g)$color <- unname(cluster_cols[tag])
  cat(sprintf("[14b] %s: %d pairs, %d nodes, %d connected components (all clonal-cluster size >=2)\n",
              tag, nrow(sub), gorder(g), n_multi))

  png(file.path(FIG_DIR, "stage_b", sprintf("14b_clonal_transmission_%s.png", tag)),
      width = 8, height = 8, units = "in", res = 300)
  par(mar = c(2, 2, 3, 2))
  set.seed(1); coords <- layout_nicely(g)
  plot(g, layout = coords, vertex.size = 8, vertex.label.cex = 0.75,
       edge.width = 2, edge.label = round(E(g)$weight, 3), edge.label.cex = 0.7,
       main = sprintf("%s clonal transmissions (fract_sites_IBD >= 0.98) — %d clonal clusters",
                      tag, n_multi))
  dev.off()
}

# also produce a summary
comp_tab <- tibble()
for (tag in unique(clonal$clonal_cluster)) {
  sub <- clonal |> filter(clonal_cluster == tag)
  g <- graph_from_data_frame(sub |> transmute(from = sample1, to = sample2), directed = FALSE)
  comp <- components(g)
  comp_tab <- bind_rows(comp_tab, tibble(
    tag = tag, n_pairs = nrow(sub), n_samples = gorder(g),
    n_clonal_clusters = sum(comp$csize > 1),
    cluster_sizes = paste(sort(comp$csize, decreasing = TRUE), collapse = ",")
  ))
}
cat("\nsummary vs poster (Mf=4, Mn=3):\n")
print(as.data.frame(comp_tab))
write_tsv(comp_tab, file.path(RES_DIR, "stage_b", "14b_clonal_cluster_counts.tsv"))
