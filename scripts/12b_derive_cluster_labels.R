# Stage B-4 helper: assign cluster labels to samples using TRUTH cleaned.3.Q
# (F-anchored per 2026-07-27), restrict to the 501 ∩ shared with regen cleaned.fam,
# and emit plink --keep + Exclude_for_{Mf,Mn,Pen}.tsv files.
suppressPackageStartupMessages({ library(tidyverse); library(here) })
source(here::here("scripts/_setup.R"))

OUT <- file.path(DATA_PROC, "cluster_maf")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Truth Q + fam (558)
truth_q  <- as.matrix(read.table(file.path(DATA_RAW, "tess3r", "cleaned.3.Q")))
truth_id <- read.table(file.path(DATA_RAW, "tess3r", "cleaned.fam"),
                       stringsAsFactors = FALSE)$V1
stopifnot(nrow(truth_q) == length(truth_id), nrow(truth_q) == 558)

# V-col → cluster from truth: V1=Mn, V2=Mf, V3=Peninsular (pinned & 100%-verified
# 2026-06-26 via Pk_clusters_metadata.csv known-membership samples).
truth_labels <- tibble(
  Sample  = truth_id,
  Cluster = c("Mn", "Mf", "Peninsular")[apply(truth_q, 1, which.max)],
  max_Q   = apply(truth_q, 1, max)
)

# Regen fam (559)
regen_id <- read.table(file.path(DATA_PROC, "plink", "cleaned.fam"),
                       stringsAsFactors = FALSE)$V1

# 501 ∩ = samples in both
inter <- intersect(regen_id, truth_id)
cat(sprintf("[12b] regen n=%d, truth n=%d, ∩=%d\n",
            length(regen_id), length(truth_id), length(inter)))

lab501 <- truth_labels |>
  filter(Sample %in% inter) |>
  arrange(match(Sample, regen_id))
cat("[12b] cluster distribution in the 501 ∩ (from truth Q):\n")
print(table(lab501$Cluster))

# Write plink --keep (FID IID = same, per --double-id convention)
write.table(data.frame(FID = lab501$Sample, IID = lab501$Sample),
            file.path(OUT, "keep_501.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Write per-sample label table
write_tsv(lab501, file.path(OUT, "cluster_labels.tsv"))

# For each cluster, Exclude_for_X.tsv = samples NOT in that cluster (so plink
# --remove keeps only that cluster). plink --remove format: FID IID (+ ok extra cols).
for (C in c("Mf", "Mn", "Peninsular")) {
  non_c <- lab501 |>
    filter(Cluster != C) |>
    transmute(FID = Sample, IID = Sample)
  out_name <- if (C == "Peninsular") "Exclude_for_Pen.tsv" else paste0("Exclude_for_", C, ".tsv")
  write.table(non_c, file.path(OUT, out_name),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  cat(sprintf("[12b] wrote %s (n=%d non-%s samples excluded)\n",
              out_name, nrow(non_c), C))
}
