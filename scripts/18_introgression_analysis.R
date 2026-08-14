# Stage B-9: introgression analysis using truth-anchored inputs.
#
# F-anchored: for validation this pass we do NOT regenerate the underlying
# hmmIBD.tsv genotype table from the no_MOI VCF (that would require running
# 17_hmmibd_genotype.R against a decompressed VCF — deferred as pipeline
# completeness). Instead we use the truth introgression outputs:
#   data/raw/Introgression/introgressed_windows_filtered.tsv   (3,943 rows)
#   data/raw/Introgression/mf_windows.tsv                      (150)
#   data/gadi/validation/introgression/mn_windows.tsv          (67; new this pass)
#   data/raw/Introgression/{mf,mn}_samples_window_counts.tsv
# and F-anchored cluster labels (data/processed/cluster_maf/cluster_labels.tsv).
#
# Reports:
#  1. Per-cluster window counts (compare to poster's "24 environment-associated").
#  2. RE-RUN the Stage-A combined-significant-windows table WITH the corrected
#     mn_windows join (the original had a copy-paste bug that read mf_windows
#     twice; now-available mn_windows fixes it). Compare to Stage-A's 34 count.
#  3. Per-sample introgression window counts, distribution.
suppressPackageStartupMessages({ library(tidyverse); library(here) })
source(here::here("scripts/_setup.R"))

INTRO_RAW <- file.path(DATA_RAW,  "Introgression")
INTRO_VAL <- file.path(PROJ, "data/gadi/validation/introgression")
OUT       <- file.path(RES_DIR, "stage_b")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Truth introgression outputs ------------------------------------------
mf_win <- read_tsv(file.path(INTRO_RAW, "mf_windows.tsv"), show_col_types = FALSE)
mn_win <- read_tsv(file.path(INTRO_VAL, "mn_windows.tsv"), show_col_types = FALSE)
filt   <- read_tsv(file.path(INTRO_RAW, "introgressed_windows_filtered.tsv"),
                   show_col_types = FALSE)
mf_samp <- read_tsv(file.path(INTRO_RAW, "mf_samples_window_counts.tsv"),
                    show_col_types = FALSE)
mn_samp <- read_tsv(file.path(INTRO_RAW, "mn_samples_window_counts.tsv"),
                    show_col_types = FALSE)

message(sprintf("[18] mf_windows: %d rows | mn_windows: %d rows | filtered: %d rows",
                nrow(mf_win), nrow(mn_win), nrow(filt)))
message(sprintf("[18] mf per-sample counts: %d rows | mn: %d rows",
                nrow(mf_samp), nrow(mn_samp)))
message(sprintf("[18] unique WINDOW ids in filtered: %d",
                length(unique(filt$WINDOW))))

# ---- Per-cluster window counts, ranges, chrom breakdown -------------------
cat("\n[18] mf_windows CHROM breakdown:\n")
print(mf_win |> count(CHROM))
cat("\n[18] mn_windows CHROM breakdown:\n")
print(mn_win |> count(CHROM))

# ---- Corrected Stage-A combined significant-windows re-run ----------------
# The original 01_pop_gen.Rmd chunk joined mf_windows.tsv twice (copy-paste
# bug in the original — should have been mf_windows + mn_windows). With
# mn_windows now available we can compute the intended union.
mf_ranges <- mf_win |>
  transmute(Window = paste0("Window_", WINDOW), CHROM, start, end)
mn_ranges <- mn_win |>
  transmute(Window = paste0("Window_", WINDOW), CHROM, start, end)
combined_ranges <- bind_rows(mf_ranges, mn_ranges) |> distinct()
message(sprintf("[18] union of mf+mn window ranges: %d rows (unique keyed by Window)",
                nrow(combined_ranges)))

# If the Stage-A significant-windows table is present, re-diff
sig_path <- file.path(RES_DIR, "sig_introgression_windows_combined.tsv")
if (file.exists(sig_path)) {
  sig <- read_tsv(sig_path, show_col_types = FALSE)
  n_sig_unique <- length(unique(sig$Window))
  n_sig_in_mf <- sum(unique(sig$Window) %in% mf_ranges$Window)
  n_sig_in_mn <- sum(unique(sig$Window) %in% mn_ranges$Window)
  message(sprintf("[18] Stage-A combined sig windows: %d unique; %d ⊂ mf_windows, %d ⊂ mn_windows",
                  n_sig_unique, n_sig_in_mf, n_sig_in_mn))
  message("[18] (poster reports 24 environment-associated windows; Stage-A found 34; ",
          "gap now decomposable into per-cluster contributions)")
  cat("[18] Stage-A sig windows per Variable (approximate — half-rows due to mf join dup):\n")
  print(sig |> distinct(Window, Variable) |> mutate(Term = sub(":.*$", "", Variable)) |>
        distinct(Window, Term) |> count(Term))
} else {
  message("[18] no Stage-A sig windows table (results/sig_introgression_windows_combined.tsv) — skip re-diff")
}

# ---- Per-sample introgression window counts vs cluster labels -------------
labels <- read_tsv(file.path(DATA_PROC, "cluster_maf", "cluster_labels.tsv"),
                   show_col_types = FALSE)
mf_samp_lab <- mf_samp |> rename(Sample = SAMPLE) |> left_join(labels, by = "Sample")
mn_samp_lab <- mn_samp |> rename(Sample = SAMPLE) |> left_join(labels, by = "Sample")
cat("\n[18] mf per-sample window counts by F-anchored cluster:\n")
print(mf_samp_lab |> count(Cluster, name = "n_samples") |> arrange(desc(n_samples)))
cat("[18] mn per-sample window counts by F-anchored cluster:\n")
print(mn_samp_lab |> count(Cluster, name = "n_samples") |> arrange(desc(n_samples)))

# ---- Write summary --------------------------------------------------------
summary_out <- tibble(
  metric = c("mf_windows_rows", "mn_windows_rows", "filtered_rows",
             "filtered_unique_windows", "union_mf_mn_windows",
             "mf_samples_with_counts", "mn_samples_with_counts"),
  value  = c(nrow(mf_win), nrow(mn_win), nrow(filt),
             length(unique(filt$WINDOW)), nrow(combined_ranges),
             nrow(mf_samp), nrow(mn_samp))
)
write_tsv(summary_out, file.path(OUT, "18_introgression_summary.tsv"))
message(sprintf("[18] wrote %s", file.path(OUT, "18_introgression_summary.tsv")))
