# Stage B-9 diagnosis: instrument the filter chain to isolate where regen
# inflates 6× vs truth (24,980 vs 3,942 rows in introgressed_windows_filtered).
# Uses the cached pre-filter dataframe (37,634 rows, 558-sample contour scan).
suppressPackageStartupMessages({ library(tidyverse); library(here) })
source(here::here("scripts/_setup.R"))

INTRO <- file.path(DATA_PROC, "introgression")
OUT   <- file.path(RES_DIR, "stage_b")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- inputs -----------------------------------------------------------------
iw_pre <- readRDS(file.path(INTRO, "introgressed_windows_prefilter.rds"))
message(sprintf("[19c] pre-filter iw: %d rows, %d unique WINDOWs, %d unique SAMPLEs, %d unique (SAMPLE,WINDOW)",
                nrow(iw_pre),
                length(unique(iw_pre$WINDOW)),
                length(unique(iw_pre$SAMPLE)),
                nrow(distinct(iw_pre, SAMPLE, WINDOW))))
cat("[19c] pre-filter cluster breakdown:\n"); print(table(iw_pre$Cluster))
cat("[19c] pre-filter (Mf/Mn/Pen) counts (Mf==0 & Mn==1) etc:\n")
print(iw_pre |> count(Cluster, Mf, Mn, Pen))

# Truth reference
truth <- read_tsv(file.path(DATA_RAW, "Introgression", "introgressed_windows_filtered.tsv"),
                  show_col_types = FALSE)
truth_windows <- unique(truth$WINDOW)
message(sprintf("[19c] TRUTH: %d rows, %d unique WINDOWs", nrow(truth), length(truth_windows)))
cat("[19c] TRUTH cluster breakdown:\n"); print(table(truth$Cluster))

# ---- step through the filter chain ----------------------------------------
report <- tibble()
step <- function(name, df, comment = "") {
  nrows <- nrow(df); nwin <- length(unique(df$WINDOW))
  message(sprintf("[19c] %-45s | rows=%6d | unique WINDOWs=%4d | %s",
                  name, nrows, nwin, comment))
  report <<- bind_rows(report, tibble(step = name, rows = nrows,
                                       unique_windows = nwin, note = comment))
}

step("0. pre-filter", iw_pre)

# STEP 1: window appears in >5 samples overall
kept <- iw_pre |> count(WINDOW) |> filter(n > 5) |> pull(WINDOW)
s1 <- iw_pre |> filter(WINDOW %in% kept)
step("1. filter WINDOW seen in >5 samples", s1)

# STEP 2: window appears in >5 samples per cluster (join-based)
per_cluster_kept <- s1 |>
  count(Cluster, WINDOW) |>
  filter(n > 5) |>
  transmute(key = paste(Cluster, WINDOW, sep = "_")) |>
  pull(key)
s2 <- s1 |>
  filter(paste(Cluster, WINDOW, sep = "_") %in% per_cluster_kept)
step("2. filter WINDOW seen in >5 samples per cluster", s2)

# STEP 3: SICAvar / KIR (from GFF) — need PKA1H1_windows to map (POS → WINDOW)
PKA1H1_windows <- read_tsv(file.path(INTRO, "PKA1H1_windows.tsv"), show_col_types = FALSE)
gff <- read.table(file.path(PROJ, "data/gadi/gff/strain_A1_H.1.Icor.gff3"),
                  sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "") |>
  rename(CHROM = V1, feature = V3, start = V4, end = V5, attribute = V9) |>
  filter(grepl("ordered", CHROM)) |>
  mutate(CHROM = str_remove(CHROM, "ordered_PKNH_"),
         CHROM = str_remove(CHROM, "_v2"))
message(sprintf("[19c] GFF rows: %d | SICA/KIR features: %d",
                nrow(gff), sum(grepl("SICA|KIR", gff$attribute))))
hyper_regions <- gff |>
  filter(grepl("SICA|KIR", attribute)) |>
  as_tibble() |>
  pivot_longer(c(start, end), names_to = "range", values_to = "POS") |>
  transmute(CHROM_POS = paste(CHROM, POS, sep = "_")) |>
  distinct() |>
  inner_join(PKA1H1_windows |> mutate(CHROM_POS = paste(CHROM, POS, sep = "_")),
             by = "CHROM_POS")
message(sprintf("[19c] SICA/KIR windows (from GFF endpoints landing in windows): %d unique WINDOW ids",
                length(unique(hyper_regions$WINDOW))))
s3 <- s2 |> filter(!(WINDOW %in% hyper_regions$WINDOW))
step("3. drop SICA/KIR windows", s3)

# STEP 4: cross-cluster hypervariable — window occurring in multiple clusters
hyper <- s3 |> distinct(WINDOW, Cluster) |> count(WINDOW) |> filter(n > 1)
s4 <- s3 |> filter(!(WINDOW %in% hyper$WINDOW))
step("4. drop windows in multiple clusters (hypervariable)", s4)

# ---- comparison vs truth --------------------------------------------------
message(sprintf("\n[19c] === POST-FILTER: regen=%d rows, %d unique WINDOWs ; TRUTH=%d rows, %d unique WINDOWs ===",
                nrow(s4), length(unique(s4$WINDOW)),
                nrow(truth), length(truth_windows)))
inter_w <- intersect(s4$WINDOW, truth_windows)
message(sprintf("[19c] WINDOW-id overlap: intersect=%d ; regen-only=%d ; truth-only=%d",
                length(inter_w),
                length(setdiff(s4$WINDOW, truth_windows)),
                length(setdiff(truth_windows, s4$WINDOW))))

# ---- experiment A: try a STRICTER "n>5" (make it "n>=5" — off-by-one?) -----
# Original R: filter(n > 5) — so windows with n=6,7,... survive. My port matches.
# BUT check: does the original perhaps sum across clusters differently?
kept_ge5 <- iw_pre |> count(WINDOW) |> filter(n >= 5) |> pull(WINDOW)
sA <- iw_pre |> filter(WINDOW %in% kept_ge5)
message(sprintf("[19c] alt-A n>=5 (off-by-one test): rows=%d ; unique WINDOWs=%d",
                nrow(sA), length(unique(sA$WINDOW))))

# ---- experiment B: what if the per-cluster >5 is intended as ONLY the win's
#                    dominant introgression-target cluster? -----------------
# Original mf_windows.tsv output uses `filter(Cluster == "Mf" & Mf == 0 & Mn == 1)`
# — i.e., "Mf sample with Mn allele". Cluster is the SAMPLE's assigned cluster.
# The pre-filter `filter(n > 5)` may have originally been AFTER an additional
# filter that also enforces the target-direction. Test: apply pre-filter within
# introgression-DIRECTIONS.

# For each sample-cluster & direction (Mf/Mn/Pen membership vector), a "direction"
# is a specific (own-cluster-allele=0, other-cluster-allele=1). Test filtering
# each direction separately then union.
dirs <- iw_pre |>
  transmute(SAMPLE, WINDOW, Cluster,
            direction = case_when(
              Cluster == "Mf" & Mf == 0 & Mn == 1 & Pen == 0 ~ "Mf<-Mn",
              Cluster == "Mn" & Mn == 0 & Mf == 1 & Pen == 0 ~ "Mn<-Mf",
              Cluster == "Peninsular" & Pen == 0 & Mf == 1 & Mn == 1 ~ "Pen<-Mf_Mn",
              TRUE ~ "OTHER"))
cat("[19c] direction breakdown of pre-filter:\n")
print(dirs |> count(direction))
# What if the count-per-window >5 was intended to be within direction?
kept_dir <- dirs |>
  filter(direction != "OTHER") |>
  count(direction, WINDOW) |>
  filter(n > 5) |>
  transmute(key = paste(direction, WINDOW, sep = "_")) |>
  pull(key)
sC <- dirs |>
  filter(paste(direction, WINDOW, sep = "_") %in% kept_dir) |>
  distinct(WINDOW, Cluster)
message(sprintf("[19c] alt-C per-direction n>5: %d unique (WINDOW,Cluster) → %d unique WINDOWs",
                nrow(sC), length(unique(sC$WINDOW))))

# ---- experiment D: what if step-4 (multi-cluster) is stricter (>= 1 = drop ANY
#                    window that ever appears in >1 cluster in the WHOLE
#                    pre-filter, not just after step 2)? --------------------
hyper_early <- iw_pre |> distinct(WINDOW, Cluster) |> count(WINDOW) |> filter(n > 1)
sD <- s3 |> filter(!(WINDOW %in% hyper_early$WINDOW))
message(sprintf("[19c] alt-D drop multi-cluster windows detected at pre-filter (not post-step-3): rows=%d ; unique WINDOWs=%d",
                nrow(sD), length(unique(sD$WINDOW))))

# ---- experiment E: apply the direction filter FIRST, then the >5 filters --
dirs2 <- iw_pre |>
  transmute(SAMPLE, WINDOW, Cluster, Mf, Mn, Pen,
            direction = case_when(
              Cluster == "Mf" & Mf == 0 & Mn == 1 & Pen == 0 ~ "Mf<-Mn",
              Cluster == "Mn" & Mn == 0 & Mf == 1 & Pen == 0 ~ "Mn<-Mf",
              Cluster == "Peninsular" & Pen == 0 & Mf == 1 & Mn == 1 ~ "Pen<-Mf_Mn",
              TRUE ~ NA_character_)) |>
  filter(!is.na(direction))
eE_step0 <- dirs2
kept_E1 <- eE_step0 |> count(WINDOW) |> filter(n > 5) |> pull(WINDOW)
sE1 <- eE_step0 |> filter(WINDOW %in% kept_E1)
kept_E2 <- sE1 |> count(direction, WINDOW) |> filter(n > 5) |>
  transmute(key = paste(direction, WINDOW, sep = "_")) |> pull(key)
sE2 <- sE1 |> filter(paste(direction, WINDOW, sep = "_") %in% kept_E2)
sE3 <- sE2 |> filter(!(WINDOW %in% hyper_regions$WINDOW))
hyper_E <- sE3 |> distinct(WINDOW, Cluster) |> count(WINDOW) |> filter(n > 1)
sE4 <- sE3 |> filter(!(WINDOW %in% hyper_E$WINDOW))
message(sprintf("[19c] alt-E direction-then-filter pipeline: rows=%d ; unique WINDOWs=%d",
                nrow(sE4), length(unique(sE4$WINDOW))))
inter_E <- intersect(sE4$WINDOW, truth_windows)
message(sprintf("[19c]         truth WINDOW overlap: %d intersect / %d truth (%d truth-only, %d regen-only)",
                length(inter_E), length(truth_windows),
                length(setdiff(truth_windows, sE4$WINDOW)),
                length(setdiff(sE4$WINDOW, truth_windows))))

write_tsv(report, file.path(OUT, "19_filter_chain_counts.tsv"))
message(sprintf("[19c] wrote %s", file.path(OUT, "19_filter_chain_counts.tsv")))
