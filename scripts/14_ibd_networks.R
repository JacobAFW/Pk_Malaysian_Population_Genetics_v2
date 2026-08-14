# Stage B-6: IBD network graphs (poster Figs 1-2 — 4 clonal clusters Mf, 3 Mn).
#
# NB: not regenerating hmmIBD locally in this pass — it takes many hours and
# requires the full VCF → hmmIBD.tsv genotype conversion; the truth fract file
# is already available and (per F decision) we anchor downstream analyses to
# poster-faithful inputs. This is documented as a caveat: hmmIBD determinism
# was NOT independently verified this pass; the network topology is truth.
#
# Reads data/raw/hmmIBD/Pk.hmm_fract_revisions.txt (11,477 pairs, 152 samples),
# joins F-anchored cluster labels from truth cleaned.3.Q, builds IBD networks
# at multiple thresholds, plots per-cluster + counts connected components as
# proxy for "clonal clusters" to compare with poster.
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(igraph)
})
source(here::here("scripts/_setup.R"))

fract <- read_tsv(file.path(DATA_RAW, "hmmIBD", "Pk.hmm_fract_revisions.txt"),
                  show_col_types = FALSE)
message(sprintf("[14] Pk.hmm_fract_revisions.txt: %d pairs, %d unique samples",
                nrow(fract), length(unique(c(fract$sample1, fract$sample2)))))

# F-anchored cluster labels (from truth cleaned.3.Q joined to cleaned.fam)
truth_q  <- as.matrix(read.table(file.path(DATA_RAW, "tess3r", "cleaned.3.Q")))
truth_id <- read.table(file.path(DATA_RAW, "tess3r", "cleaned.fam"),
                       stringsAsFactors = FALSE)$V1
truth_labels <- tibble(
  Sample  = truth_id,
  Cluster = c("Mn", "Mf", "Peninsular")[apply(truth_q, 1, which.max)]
)

# Sample-side cluster coverage in the hmm fract file
samples_in_fract <- unique(c(fract$sample1, fract$sample2))
lab_frac <- truth_labels |> filter(Sample %in% samples_in_fract)
message(sprintf("[14] fract-file samples with F-anchored labels: %d / %d",
                nrow(lab_frac), length(samples_in_fract)))
cat("[14] cluster distribution in fract-file samples (labelled subset):\n")
print(table(lab_frac$Cluster))

# ---- IBD threshold + network builder ---------------------------------------
build_graph <- function(fract_df, thr) {
  edges <- fract_df |>
    filter(fract_sites_IBD >= thr) |>
    transmute(from = sample1, to = sample2, weight = fract_sites_IBD)
  g <- igraph::graph_from_data_frame(edges,
                                     vertices = data.frame(name = samples_in_fract),
                                     directed = FALSE)
  # attach cluster attribute
  vc <- truth_labels$Cluster[match(V(g)$name, truth_labels$Sample)]
  vc[is.na(vc)] <- "Unlabelled"
  V(g)$Cluster <- vc
  V(g)$color   <- unname(c(cluster_cols, "Unlabelled" = "grey60")[vc])
  g
}

# Poster-relevant cutoffs — 0.20 (broad "significant IBD") + 0.50 (clonal)
cutoffs <- c(0.20, 0.50)
cc_report <- tibble()
for (thr in cutoffs) {
  g <- build_graph(fract, thr)
  # Full network plot (all clusters)
  set.seed(1); coords <- layout_nicely(g)
  png(file.path(FIG_DIR, "stage_b", sprintf("14_ibd_network_all_thr%.2f.png", thr)),
      width = 10, height = 10, units = "in", res = 300)
  par(mar = c(2, 2, 3, 2))
  plot(g, layout = coords, vertex.size = 3, vertex.label = NA,
       edge.width = 0.5, edge.color = "grey70",
       main = sprintf("IBD network (all clusters, fract_sites_IBD >= %.2f)", thr))
  legend("bottomright", legend = names(cluster_cols),
         fill = unname(cluster_cols), bty = "n", cex = 0.9)
  dev.off()

  # Per-cluster subgraph + connected components
  for (C in c("Mf", "Mn", "Peninsular")) {
    keep <- V(g)$Cluster == C
    sub_g <- induced_subgraph(g, which(keep))
    # only edges within the cluster
    comps <- components(sub_g)
    n_comp   <- comps$no
    n_multi  <- sum(comps$csize > 1)                      # non-singleton components
    n_vert   <- gorder(sub_g)
    n_edges  <- gsize(sub_g)
    cc_report <- bind_rows(cc_report, tibble(
      threshold = thr, cluster = C,
      n_samples = n_vert, n_edges = n_edges,
      n_components = n_comp, n_clonal_clusters = n_multi
    ))
    png(file.path(FIG_DIR, "stage_b", sprintf("14_ibd_network_%s_thr%.2f.png", C, thr)),
        width = 8, height = 8, units = "in", res = 300)
    par(mar = c(2, 2, 3, 2))
    set.seed(1); coords_c <- layout_nicely(sub_g)
    plot(sub_g, layout = coords_c, vertex.size = 4, vertex.label = NA,
         edge.width = 0.6, edge.color = "grey70",
         vertex.color = unname(cluster_cols[C]),
         main = sprintf("%s IBD network (fract_sites_IBD >= %.2f) — %d nodes / %d edges / %d clonal clusters",
                        C, thr, n_vert, n_edges, n_multi))
    dev.off()
  }
}
cat("\n[14] clonal-cluster counts (connected components with >1 vertex):\n")
print(as.data.frame(cc_report))
write_tsv(cc_report, file.path(RES_DIR, "stage_b", "14_ibd_component_counts.tsv"))
message(sprintf("[14] wrote %s and %d figures under figures/stage_b/",
                file.path(RES_DIR, "stage_b", "14_ibd_component_counts.tsv"),
                length(cutoffs) * 4))
