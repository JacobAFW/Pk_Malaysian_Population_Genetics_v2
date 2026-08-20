# Stage B-8 v2 helper: derive REGEN ADMIXTURE cluster labels for the full 559
# cohort (drop the 501 ∩ F-anchor for stage 8 — selection needs cohort SIZE, not
# poster-label purity). V-col → cluster identity pinned from known-membership
# metadata (Pk_clusters_metadata.csv + Pk_clusters_peninsular_metadata.csv);
# same technique as the Stage-A Fix-1/2 helper.
suppressPackageStartupMessages({ library(tidyverse); library(here) })
source(here::here("scripts/_setup.R"))

OUT <- file.path(DATA_PROC, "selection_v2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# regen ADMIXTURE Q + fam
q  <- as.matrix(read.table(file.path(DATA_PROC, "admixture", "cleaned.3.Q")))
id <- read.table(file.path(DATA_PROC, "plink", "cleaned.fam"),
                 stringsAsFactors = FALSE)$V1
stopifnot(nrow(q) == length(id), nrow(q) == 559)

# known-membership metadata
km <- read.csv(file.path(PROJ, "data/gadi/Pk_clusters_metadata.csv"),
               strip.white = TRUE)
km <- data.frame(Sample = trimws(km$Sample),
                 Cluster = sub("-Pk$", "", trimws(km$Group)))
kp <- read.csv(file.path(PROJ, "data/gadi/Pk_clusters_peninsular_metadata.csv"),
               strip.white = TRUE)
known <- rbind(km,
               data.frame(Sample = trimws(kp$ENA_accession_no_ER),
                          Cluster = "Peninsular"))
known <- known[known$Cluster %in% c("Peninsular","Mf","Mn"), ]

# V-col → cluster from known argmax
qdf <- data.frame(Sample = id, V1 = q[,1], V2 = q[,2], V3 = q[,3])
qdf_lab <- merge(qdf, known, by = "Sample")
argmax <- apply(qdf_lab[, c("V1","V2","V3")], 1, which.max)
tab <- table(Vcol = factor(paste0("V", argmax), levels = c("V1","V2","V3")),
             Cluster = qdf_lab$Cluster)
cat("[16v2] Q-col vs known cluster tally (n=", nrow(qdf_lab), " known):\n", sep = "")
print(tab)

col_to_cluster <- character(3); names(col_to_cluster) <- c("V1","V2","V3")
row_sums <- rowSums(tab)
for (v in c("V1","V2","V3")) {
  if (row_sums[[v]] > 0) col_to_cluster[v] <- names(which.max(tab[v, ]))
}
empty_v <- names(col_to_cluster)[col_to_cluster == ""]
if (length(empty_v) == 1) {
  missing_c <- setdiff(c("Peninsular","Mf","Mn"), col_to_cluster)
  if (length(missing_c) == 1) col_to_cluster[empty_v] <- missing_c
}
cat("[16v2] V-col → cluster:", paste(names(col_to_cluster), col_to_cluster, sep = "=", collapse = ", "), "\n")

# assign regen labels to ALL 559 samples via argmax
labels <- tibble(
  Sample  = id,
  Cluster = col_to_cluster[apply(q, 1, which.max)],
  max_Q   = apply(q, 1, max)
)
cat("[16v2] regen 559 cluster distribution (from ADMIXTURE argmax):\n")
print(table(labels$Cluster))

# Write labels + per-cluster exclude files.
# Two formats: plink --remove wants FID\tIID (2 col); bcftools -S wants 1 col.
write_tsv(labels, file.path(OUT, "regen_labels.tsv"))
for (C in c("Mf", "Mn", "Peninsular")) {
  non_c <- labels |> filter(Cluster != C) |> pull(Sample)
  tag <- if (C == "Peninsular") "Pen" else C

  # plink format (2-col)
  write.table(data.frame(FID = non_c, IID = non_c),
              file.path(OUT, paste0("Exclude_for_", tag, ".tsv")),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  # bcftools format (1-col, one sample per line)
  writeLines(non_c, file.path(OUT, paste0("Exclude_for_", tag, "_bcftools.txt")))

  cat(sprintf("[16v2] wrote Exclude_for_%s.tsv (2-col plink) + _bcftools.txt (1-col): %d non-%s samples\n",
              tag, length(non_c), C))
}
