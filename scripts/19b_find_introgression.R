# Stage B-9: full introgression pipeline — faithful lift of
# scripts/gadi/find_introgression_updated_framework.R with local paths.
# F-anchored: uses truth cluster labels (from data/raw/tess3r/cleaned.3.Q).
#
# NB (per prompt 2026-08-17): provisional — Jacob is rewriting this script to be
# data-agnostic, so this pass is a faithful lift, don't over-invest in matching.
# Regenerate mn_windows (was absent locally at Stage-A), diff vs truth files,
# report 34-vs-24 with the no-FDR caveat, don't "fix" it.
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(data.table); library(sp); library(MASS)
})
source(here::here("scripts/_setup.R"))

INTRO <- file.path(DATA_PROC, "introgression")
dir.create(INTRO, showWarnings = FALSE, recursive = TRUE)

# ---- inputs ---------------------------------------------------------------
# F-anchored labels — truth 558 cluster assignments, joined against the samples
# present in the hmmIBD.tsv genotype table (whole regen cohort).
truth_q  <- as.matrix(read.table(file.path(DATA_RAW, "tess3r", "cleaned.3.Q")))
truth_id <- read.table(file.path(DATA_RAW, "tess3r", "cleaned.fam"),
                       stringsAsFactors = FALSE)$V1
metadata <- tibble(Sample  = truth_id,
                   Cluster = c("Mn","Mf","Peninsular")[apply(truth_q, 1, which.max)])

genotype_table <- read_tsv(file.path(INTRO, "hmmIBD.tsv"),
                           show_col_types = FALSE, guess_max = 1e5)
message(sprintf("[19b] hmmIBD.tsv: %d rows × %d cols  (cols 3.. are samples)",
                nrow(genotype_table), ncol(genotype_table)))

# Only keep sample columns whose ID is in truth metadata (drops the ~200 samples
# with no truth cluster label — F decision).
sample_cols <- setdiff(names(genotype_table), c("CHROM", "POS"))
keep_cols <- sample_cols[sample_cols %in% metadata$Sample]
message(sprintf("[19b] samples in hmmIBD.tsv: %d ; kept (F-anchored): %d",
                length(sample_cols), length(keep_cols)))
genotype_table <- genotype_table[, c("CHROM","POS", keep_cols)]

# Long form
genotype_table_long <- genotype_table |>
  pivot_longer(3:ncol(genotype_table), names_to = "SAMPLE", values_to = "SNP") |>
  mutate(SNP = as.numeric(SNP)) |>
  left_join(metadata |> rename(SAMPLE = Sample), by = "SAMPLE")

# Dominant allele per cluster per site
dominant_allele <- genotype_table_long |>
  filter(SNP >= 0) |>
  mutate(SNP = as.factor(SNP)) |>
  group_by(CHROM, POS, Cluster, SNP) |>
  summarise(Allele_count = n(), .groups = "drop") |>
  group_by(CHROM, POS, Cluster) |>
  filter(Allele_count == max(Allele_count)) |>
  ungroup() |>
  dplyr::select(-Allele_count) |>
  group_by(CHROM, POS) |>
  mutate(row = cur_group_id()) |>
  drop_na() |>
  group_by(row) |>
  # Old dplyr recycled here; modern requires reframe(). Emits one row per
  # (CHROM,POS,Cluster,SNP) combo within each row-group, with n = group size.
  reframe(CHROM = CHROM, POS = POS, Cluster = Cluster, SNP = SNP, n = n()) |>
  filter(n < 4) |>
  pivot_wider(names_from = Cluster, values_from = SNP) |>
  dplyr::select(-c(row, n))

message(sprintf("[19b] dominant_allele: %d sites", nrow(dominant_allele)))

# introgression_table (Intro_clust label per sample × site) — nested-ifelse
# ports of the original (see find_introgression_updated_framework.R L48+).
introgression_table <- genotype_table_long |>
  left_join(dominant_allele, by = c("CHROM","POS")) |>
  mutate(SNP = as.character(SNP)) |>
  mutate(across(c(Mn, Mf, Peninsular), as.character)) |>
  mutate(Intro_clust = SNP) |>
  mutate(Intro_clust = case_when(
    Cluster == "Mn" & Intro_clust == Mn ~ "Mn",
    Cluster == "Mn" & Intro_clust == Mf  & Intro_clust != Peninsular ~ "Mf",
    Cluster == "Mn" & Intro_clust == Peninsular & Intro_clust != Mf ~ "Peninsular",
    Cluster == "Mn" & Intro_clust == Mf  & Intro_clust == Peninsular ~ "Mf_Peninsular",
    Cluster == "Mn" ~ "Mn",
    Cluster == "Mf" & Intro_clust == Mf ~ "Mf",
    Cluster == "Mf" & Intro_clust == Mn  & Intro_clust != Peninsular ~ "Mn",
    Cluster == "Mf" & Intro_clust == Peninsular & Intro_clust != Mn ~ "Peninsular",
    Cluster == "Mf" & Intro_clust == Mn  & Intro_clust == Peninsular ~ "Mn_Peninsular",
    Cluster == "Mf" ~ "Mf",
    Cluster == "Peninsular" & Intro_clust == Peninsular ~ "Peninsular",
    Cluster == "Peninsular" & Intro_clust == Mf & Intro_clust != Mn ~ "Mf",
    Cluster == "Peninsular" & Intro_clust == Mn & Intro_clust != Mf ~ "Mn",
    Cluster == "Peninsular" & Intro_clust == Mf & Intro_clust == Mn ~ "Mf_Mn",
    Cluster == "Peninsular" ~ "Peninsular",
    TRUE ~ Intro_clust
  ))

write_tsv(introgression_table, file.path(INTRO, "introgression_table.tsv"))
message(sprintf("[19b] introgression_table: %d rows", nrow(introgression_table)))

# 10kb sliding-window genetic distances
window_size <- 10000L
bed <- read.table(file.path(PROJ, "data/gadi/bed/strain_A1_H.1.Icor.fasta.bed"),
                  stringsAsFactors = FALSE)
PKA1H1_windows <- bed |>
  as_tibble() |>
  rename(CHROM = V1, Start = V2, End = V3) |>
  mutate(CHROM = str_remove(CHROM, "ordered_PKNH_"),
         CHROM = str_remove(CHROM, "_v2")) |>
  filter(CHROM != "PKNH_MIT", CHROM != "new_API_strain_A1_H.1")

# expand to per-POS windows
suppressPackageStartupMessages(library(plyr))
PKA1H1_windows <- plyr::ddply(PKA1H1_windows, "CHROM",
                              summarise, POS = seq(Start, End)) |>
  as_tibble() |>
  mutate(TMP_WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) |>
  # Old dplyr: group_indices(group_by(...)); modern: dense_rank on the key.
  mutate(WINDOW = dplyr::dense_rank(paste(CHROM, TMP_WINDOW, sep = "_"))) |>
  dplyr::select(-TMP_WINDOW)
detach("package:plyr", unload = TRUE)
write_tsv(PKA1H1_windows, file.path(INTRO, "PKA1H1_windows.tsv"))
message(sprintf("[19b] PKA1H1_windows: %d POS rows, %d unique windows",
                nrow(PKA1H1_windows), length(unique(PKA1H1_windows$WINDOW))))

# sample × window distance table
introgression_table_window <- introgression_table |>
  mutate(POS = as.numeric(POS)) |>
  left_join(PKA1H1_windows, by = c("CHROM","POS")) |>
  filter(SNP >= 0, !is.na(WINDOW)) |>
  group_by(SAMPLE, WINDOW) |>
  mutate(n = n()) |>
  filter(n > 5) |>
  summarise(
    Mf_distance  = sum(SNP != Mf,          na.rm = TRUE) / n() * 100,
    Mn_distance  = sum(SNP != Mn,          na.rm = TRUE) / n() * 100,
    Pen_distance = sum(SNP != Peninsular,  na.rm = TRUE) / n() * 100,
    .groups = "drop"
  ) |>
  left_join(metadata |> rename(SAMPLE = Sample), by = "SAMPLE") |>
  filter(SAMPLE %in% metadata$Sample) |>
  filter(Mn_distance != Mf_distance)

write_tsv(introgression_table_window, file.path(INTRO, "introgression_table_window.tsv"))
message(sprintf("[19b] introgression_table_window: %d rows",
                nrow(introgression_table_window)))

# ---- contour-based introgression detection --------------------------------
find_introgressed_regions <- function(SAMPLENAME, itw = introgression_table_window) {
  raster_plot <- itw |>
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = itw |> filter(SAMPLE == SAMPLENAME), size = 0.75) +
    geom_density_2d(mapping = aes(x = Mf_distance, y = Mn_distance, colour = Cluster),
                    data = itw, contour_var = "density") +
    scale_colour_manual(values = c("Mf" = "#440154FF", "Mn" = "#73D055FF",
                                    "Peninsular" = "#39568CFF")) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100))

  pb <- ggplot_build(raster_plot)
  pts <- pb$data[[1]]
  contours <- pb$data[[2]]

  make_membership <- function(colour_hex) {
    c <- contours |> filter(colour == colour_hex, level > 5e-4)
    if (nrow(c) == 0) return(rep(0, nrow(pts)))
    # Original used positional c[[3]]/c[[4]] which in modern ggplot2 (>=3.4)
    # maps to `y` and `piece` — nonsense as polygon coords. Use named cols.
    point.in.polygon(pol.x = c$x, pol.y = c$y,
                     point.x = pts$x, point.y = pts$y)
  }

  Mf_in  <- make_membership("#440154FF")
  Mn_in  <- make_membership("#73D055FF")
  Pen_in <- make_membership("#39568CFF")

  itw |>
    filter(SAMPLE == SAMPLENAME) |>
    mutate(Mf = Mf_in, Mn = Mn_in, Pen = Pen_in) |>
    filter(case_when(
      Cluster == "Mn"          ~ Mf == 1 & Mn == 0 & Pen == 0,
      Cluster == "Mf"          ~ Mf == 0 & Mn == 1 & Pen == 0,
      Cluster == "Peninsular"  ~ Mf == 1 & Mn == 1 & Pen == 0,
      TRUE ~ FALSE))
}

sample_names <- unique(introgression_table_window$SAMPLE)
.iw_cache <- file.path(INTRO, "introgressed_windows_prefilter.rds")
if (file.exists(.iw_cache)) {
  message(sprintf("[19b] loading cached pre-filter introgressed_windows from %s", .iw_cache))
  introgressed_windows <- readRDS(.iw_cache)
} else {
  message(sprintf("[19b] scanning %d samples for introgressed windows...", length(sample_names)))
  introgressed_windows <- purrr::map_dfr(sample_names, function(s) {
    tryCatch(find_introgressed_regions(s), error = function(e) tibble())
  })
  saveRDS(introgressed_windows, .iw_cache)
}
message(sprintf("[19b] introgressed_windows (pre-filter): %d rows", nrow(introgressed_windows)))

# ---- filter chain (as in original) ----------------------------------------
# n>5 per window overall
kept <- introgressed_windows |> count(WINDOW) |> filter(n > 5) |> pull(WINDOW)
introgressed_windows <- introgressed_windows |> filter(WINDOW %in% kept)

# n>5 per window per cluster
per_cluster_kept <- introgressed_windows |>
  count(Cluster, WINDOW) |>
  filter(n > 5) |>
  transmute(key = paste(Cluster, WINDOW, sep = "_")) |>
  pull(key)
introgressed_windows <- introgressed_windows |>
  filter(paste(Cluster, WINDOW, sep = "_") %in% per_cluster_kept)

# SICAvar / KIR filter
# GFF has ## header lines + ### record separators; leave comment.char="#" default
# so read.table skips both. Use fill=TRUE for robustness.
gff <- read.table(file.path(PROJ, "data/gadi/gff/strain_A1_H.1.Icor.gff3"),
                  sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "") |>
  rename(CHROM = V1, feature = V3, start = V4, end = V5, attribute = V9) |>
  filter(grepl("ordered", CHROM)) |>
  mutate(CHROM = str_remove(CHROM, "ordered_PKNH_"),
         CHROM = str_remove(CHROM, "_v2"))
hyper_regions <- gff |>
  filter(grepl("SICA|KIR", attribute)) |>
  as_tibble() |>
  pivot_longer(c(start, end), names_to = "range", values_to = "POS") |>
  transmute(CHROM_POS = paste(CHROM, POS, sep = "_")) |>
  distinct() |>
  inner_join(PKA1H1_windows |> mutate(CHROM_POS = paste(CHROM, POS, sep = "_")),
             by = "CHROM_POS")
introgressed_windows <- introgressed_windows |>
  filter(!(WINDOW %in% hyper_regions$WINDOW))

# hypervariable (windows in multiple clusters)
hyper <- introgressed_windows |> distinct(WINDOW, Cluster) |> count(WINDOW) |> filter(n > 1)
introgressed_windows <- introgressed_windows |> filter(!(WINDOW %in% hyper$WINDOW))

write_tsv(introgressed_windows, file.path(INTRO, "introgressed_windows_filtered.tsv"))
message(sprintf("[19b] introgressed_windows_filtered: %d rows", nrow(introgressed_windows)))

# ---- outputs matching truth format ----------------------------------------
window_ranges <- PKA1H1_windows |>
  group_by(CHROM, WINDOW) |>
  summarise(start = min(POS), end = max(POS), .groups = "drop")

mf_windows <- introgressed_windows |>
  filter(Cluster == "Mf" & Mf == 0 & Mn == 1) |>
  count(WINDOW) |>
  left_join(window_ranges, by = "WINDOW") |>
  arrange(desc(n))
write_tsv(mf_windows, file.path(INTRO, "mf_windows.tsv"))

mn_windows <- introgressed_windows |>
  filter(Cluster == "Mn" & Mf == 1 & Mn == 0) |>
  count(WINDOW) |>
  left_join(window_ranges, by = "WINDOW") |>
  arrange(desc(n))
write_tsv(mn_windows, file.path(INTRO, "mn_windows.tsv"))

mf_samples <- introgressed_windows |>
  filter(Cluster == "Mf" & Mf == 0 & Mn == 1) |>
  count(SAMPLE) |> arrange(desc(n))
write_tsv(mf_samples, file.path(INTRO, "mf_samples_window_counts.tsv"))

mn_samples <- introgressed_windows |>
  filter(Cluster == "Mn" & Mf == 1 & Mn == 0) |>
  count(SAMPLE) |> arrange(desc(n))
write_tsv(mn_samples, file.path(INTRO, "mn_samples_window_counts.tsv"))

message(sprintf("[19b] outputs: mf_windows=%d rows, mn_windows=%d rows, mf_samples=%d, mn_samples=%d",
                nrow(mf_windows), nrow(mn_windows), nrow(mf_samples), nrow(mn_samples)))

# ---- diff vs truth --------------------------------------------------------
diff_rows <- tibble(
  file = c("mf_windows.tsv", "mn_windows.tsv",
           "mf_samples_window_counts.tsv", "mn_samples_window_counts.tsv",
           "introgressed_windows_filtered.tsv"),
  regen_n = c(nrow(mf_windows), nrow(mn_windows),
              nrow(mf_samples), nrow(mn_samples), nrow(introgressed_windows)),
  truth_source = c("data/raw/Introgression/mf_windows.tsv",
                   "data/gadi/validation/introgression/mn_windows.tsv",
                   "data/raw/Introgression/mf_samples_window_counts.tsv",
                   "data/raw/Introgression/mn_samples_window_counts.tsv",
                   "data/raw/Introgression/introgressed_windows_filtered.tsv"),
  truth_n = c(nrow(read_tsv(file.path(DATA_RAW, "Introgression", "mf_windows.tsv"), show_col_types=FALSE)),
              nrow(read_tsv(file.path(PROJ, "data/gadi/validation/introgression", "mn_windows.tsv"), show_col_types=FALSE)),
              nrow(read_tsv(file.path(DATA_RAW, "Introgression", "mf_samples_window_counts.tsv"), show_col_types=FALSE)),
              nrow(read_tsv(file.path(DATA_RAW, "Introgression", "mn_samples_window_counts.tsv"), show_col_types=FALSE)),
              nrow(read_tsv(file.path(DATA_RAW, "Introgression", "introgressed_windows_filtered.tsv"), show_col_types=FALSE)))
)
diff_rows$diff <- diff_rows$regen_n - diff_rows$truth_n
cat("\n[19b] regen vs truth diff:\n"); print(as.data.frame(diff_rows))
write_tsv(diff_rows, file.path(RES_DIR, "stage_b", "19_introgression_diff.tsv"))
message(sprintf("[19b] wrote %s", file.path(RES_DIR, "stage_b", "19_introgression_diff.tsv")))
