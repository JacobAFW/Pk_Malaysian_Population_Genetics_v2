# Stage B-3 follow-up: argmax-agreement vs truth restricted to CLUSTER CORES
# (samples where truth max-Q > threshold) — to test whether the 88% disagreement
# is concentrated at the genuinely admixed boundary (poster: "strong Mf-Mn
# intermixing") rather than at the cluster cores.
suppressPackageStartupMessages({ library(tidyverse); library(here); library(gtools) })
source(here::here("scripts/_setup.R"))

regen_q  <- as.matrix(read.table(file.path(DATA_PROC, "admixture", "cleaned.3.Q")))
regen_id <- read.table(file.path(DATA_PROC, "plink", "cleaned.fam"), stringsAsFactors = FALSE)$V1
truth_q  <- as.matrix(read.table(file.path(DATA_RAW,  "tess3r",    "cleaned.3.Q")))
truth_id <- read.table(file.path(DATA_RAW,  "tess3r",    "cleaned.fam"), stringsAsFactors = FALSE)$V1
stopifnot(ncol(regen_q) == 3, ncol(truth_q) == 3)

inter <- intersect(regen_id, truth_id)
r <- regen_q[match(inter, regen_id), , drop = FALSE]
t <- truth_q [match(inter, truth_id), , drop = FALSE]

# pick best permutation by mean per-sample correlation (same as 11b)
perms <- gtools::permutations(3, 3)
best_p <- NULL; best_mc <- -Inf
for (i in seq_len(nrow(perms))) {
  p <- perms[i, ]
  mc <- mean(vapply(seq_len(nrow(r)), function(s) cor(r[s, p], t[s, ]), numeric(1)), na.rm = TRUE)
  if (mc > best_mc) { best_p <- p; best_mc <- mc }
}
r <- r[, best_p]

argmax_r <- apply(r, 1, which.max)
argmax_t <- apply(t, 1, which.max)
maxq_t   <- apply(t, 1, max)

# Helper: agreement at varying truth-Q thresholds
cat("argmax agreement at truth-max-Q thresholds (501 ∩ samples):\n\n")
cat(sprintf("%-12s %-10s %-12s %s\n", "threshold", "n_core", "n_agree", "agreement"))
for (thr in c(0.50, 0.70, 0.80, 0.90, 0.95, 0.99)) {
  mask <- maxq_t >= thr
  n_core <- sum(mask)
  n_agree <- sum(argmax_r[mask] == argmax_t[mask])
  cat(sprintf("%-12.2f %-10d %-12d %.4f\n", thr, n_core, n_agree, n_agree / n_core))
}

cat("\nbreakdown at threshold 0.95 (CORE samples per cluster, truth-side argmax):\n")
mask95 <- maxq_t >= 0.95
core_truth_cluster <- argmax_t[mask95]
core_regen_cluster <- argmax_r[mask95]
print(table(truth = core_truth_cluster, regen = core_regen_cluster))

# Re-map V-col → cluster names using the same metadata pin as 01_pop_gen.Rmd
known_mfmn <- read.csv(file.path(PROJ, "data/gadi/Pk_clusters_metadata.csv"),
                       strip.white = TRUE)
known_mfmn <- data.frame(Sample  = trimws(known_mfmn$Sample),
                         Cluster = sub("-Pk$", "", trimws(known_mfmn$Group)))
known_pen <- read.csv(file.path(PROJ, "data/gadi/Pk_clusters_peninsular_metadata.csv"),
                      strip.white = TRUE)
known <- rbind(known_mfmn,
               data.frame(Sample = trimws(known_pen$ENA_accession_no_ER),
                          Cluster = "Peninsular"))
known <- known[known$Cluster %in% c("Peninsular","Mf","Mn"), ]

# Truth Q with sample IDs
truth_named <- data.frame(Sample = truth_id, V1 = truth_q[,1], V2 = truth_q[,2], V3 = truth_q[,3])
join_t <- merge(truth_named, known, by = "Sample")
cat(sprintf("\nknown-cluster samples in truth Q (for V-col→name pin): %d\n", nrow(join_t)))
arg_t_known <- apply(join_t[, c("V1","V2","V3")], 1, which.max)
print(table(V = paste0("V", arg_t_known), Cluster = join_t$Cluster))
