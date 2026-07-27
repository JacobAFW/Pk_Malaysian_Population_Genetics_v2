# Stage B-2 validator: diff regenerated cleaned.fam vs data/raw/tess3r/cleaned.fam
# (sample set + sample IDs only; metadata columns differ by design).
suppressPackageStartupMessages({ library(tidyverse); library(here) })
source(here::here("scripts/_setup.R"))

regen <- read.table(file.path(DATA_PROC, "plink",  "cleaned.fam"),
                    stringsAsFactors = FALSE)$V1
truth <- read.table(file.path(DATA_RAW,  "tess3r", "cleaned.fam"),
                    stringsAsFactors = FALSE)$V1

cat(sprintf("regen n = %d   |   truth n = %d   (poster target = 558)\n",
            length(regen), length(truth)))

inter   <- intersect(regen, truth)
only_r  <- setdiff(regen, truth)
only_t  <- setdiff(truth, regen)
cat(sprintf("intersect = %d   |   regen-only = %d   |   truth-only = %d\n",
            length(inter), length(only_r), length(only_t)))

rep <- file.path(RES_DIR, "stage_b", "10_cleaned_fam_diff.tsv")
dir.create(dirname(rep), showWarnings = FALSE, recursive = TRUE)
tibble(Sample = c(only_r, only_t),
       Side   = c(rep("regen_only", length(only_r)),
                  rep("truth_only", length(only_t)))) |>
  write_tsv(rep)
cat(sprintf("wrote %s (n=%d divergent samples)\n", rep, length(only_r) + length(only_t)))

if (length(regen) != 558) {
  cat(sprintf("\n*** CHECKPOINT FAIL: cleaned.fam = %d, not 558 ***\n", length(regen)))
} else if (length(only_r) == 0 && length(only_t) == 0) {
  cat("\n*** CHECKPOINT PASS: 558 samples, sample set matches truth exactly ***\n")
} else {
  cat(sprintf("\n*** CHECKPOINT PARTIAL: n=558 but sample set differs (%d divergent) — see %s ***\n",
              length(only_r) + length(only_t), rep))
}
