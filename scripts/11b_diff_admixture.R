# Stage B-3 validator: diff regenerated cleaned.3.Q vs data/raw/tess3r/cleaned.3.Q
# by per-sample cluster assignment (permutation-invariant), since ADMIXTURE
# label-switches columns across runs.
suppressPackageStartupMessages({ library(tidyverse); library(here); library(gtools) })
source(here::here("scripts/_setup.R"))

regen_q  <- as.matrix(read.table(file.path(DATA_PROC, "admixture", "cleaned.3.Q")))
regen_id <- read.table(file.path(DATA_PROC, "plink", "cleaned.fam"), stringsAsFactors = FALSE)$V1

truth_q  <- as.matrix(read.table(file.path(DATA_RAW,  "tess3r",    "cleaned.3.Q")))
truth_id <- read.table(file.path(DATA_RAW,  "tess3r",    "cleaned.fam"), stringsAsFactors = FALSE)$V1

stopifnot(ncol(regen_q) == 3, ncol(truth_q) == 3)
cat(sprintf("regen: n=%d (truth: n=%d)\n", length(regen_id), length(truth_id)))

# Sample-set overlap
inter <- intersect(regen_id, truth_id)
cat(sprintf("sample intersection = %d  |  regen-only = %d  |  truth-only = %d\n",
            length(inter), length(setdiff(regen_id, truth_id)), length(setdiff(truth_id, regen_id))))

if (length(inter) == 0) stop("no overlapping samples — cannot diff")

# Align by Sample
r <- regen_q[match(inter, regen_id), , drop = FALSE]
t <- truth_q [match(inter, truth_id), , drop = FALSE]

# Try all 6 column permutations; pick the one with max mean per-sample correlation
perms <- gtools::permutations(3, 3)
best <- list(perm = NULL, mean_cor = -Inf, conf = NULL, agree = NULL)
for (i in seq_len(nrow(perms))) {
  p <- perms[i, ]
  r_p <- r[, p]
  per_sample_cor <- vapply(seq_len(nrow(r_p)),
                           function(s) cor(r_p[s, ], t[s, ]),
                           numeric(1))
  argmax_r <- apply(r_p, 1, which.max)
  argmax_t <- apply(t,   1, which.max)
  agree <- mean(argmax_r == argmax_t)
  mc <- mean(per_sample_cor, na.rm = TRUE)
  if (mc > best$mean_cor) {
    best <- list(perm = p, mean_cor = mc, agree = agree,
                 argmax_r = argmax_r, argmax_t = argmax_t,
                 r_aligned = r_p)
  }
}
cat(sprintf("best column permutation (regen col -> truth col): %s\n",
            paste(seq_len(3), "->", best$perm, collapse = ", ")))
cat(sprintf("mean per-sample Q correlation (after permutation): %.4f\n", best$mean_cor))
cat(sprintf("argmax-cluster agreement (after permutation): %.4f  (%d / %d samples)\n",
            best$agree, sum(best$argmax_r == best$argmax_t), nrow(r)))

cat("\nconfusion table (rows = regen col after permutation, cols = truth argmax col):\n")
print(table(regen = best$argmax_r, truth = best$argmax_t))

# Per-sample absolute Q difference (after permutation)
diffmat <- abs(best$r_aligned - t)
cat(sprintf("\nper-cell |Q_regen - Q_truth| summary: mean=%.4f  median=%.4f  max=%.4f\n",
            mean(diffmat), median(diffmat), max(diffmat)))
cat(sprintf("samples with max-Q diff > 0.05 across cells: %d / %d\n",
            sum(apply(diffmat, 1, max) > 0.05), nrow(diffmat)))

# Write a small report
rep_path <- file.path(RES_DIR, "stage_b", "11_admixture_diff.tsv")
dir.create(dirname(rep_path), showWarnings = FALSE, recursive = TRUE)
tibble(Sample = inter,
       argmax_regen = best$argmax_r,
       argmax_truth = best$argmax_t,
       agree = (best$argmax_r == best$argmax_t),
       max_abs_diff = apply(diffmat, 1, max),
       sum_abs_diff = rowSums(diffmat)) |>
  write_tsv(rep_path)
cat(sprintf("wrote %s\n", rep_path))
