# Pipeline-completeness (not executed for validation this pass): VCF → hmmIBD.tsv.
# Ports scripts/gadi/hmmIBD_genotype_file.R with local paths + parametric inputs.
#
# hmmIBD input format: tab-delimited, one SNP per row:
#   chrom  pos  sample1  sample2  ...
# Genotype codes: 0 = ref allele, 1..4 = alt alleles, -1 = missing.
#
# Reads:  VCF (bcftools view -R snp_positions -S ^exclude ...) → gt/GT
# Writes: hmmIBD.tsv (whole-cohort) or Mf/Mn/Pen_hmmIBD.tsv (cluster-spec).
#
# Usage (R):
#   source("scripts/17_hmmibd_genotype.R")
#   make_hmmibd_tsv(vcf_path, grep_patterns_path, out_tsv)
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })

make_hmmibd_tsv <- function(vcf_path, grep_patterns_path, out_tsv, skip_header_lines = 72) {
  stopifnot(file.exists(vcf_path), file.exists(grep_patterns_path))
  # skip_header_lines: number of ## header lines in the VCF before #CHROM row.
  # (72 in the original, but varies — see wrapper below that autodetects.)
  read_tsv(grep_patterns_path, col_names = c("CHROM", "POS"), show_col_types = FALSE) |>
    mutate(POS = as.numeric(POS)) |>
    left_join(
      read_table(vcf_path, skip = skip_header_lines, show_col_types = FALSE) |>
        rename(CHROM = `#CHROM`) |> mutate(POS = as.numeric(POS)),
      by = c("CHROM", "POS")
    ) |>
    mutate(CHROM = str_remove(CHROM, "ordered_PKNH_"),
           CHROM = str_remove(CHROM, "_v2")) |>
    mutate_at(c(10:ncol(across(everything()))), ~ str_remove(., ":.*")) |>
    mutate_at(c(10:ncol(across(everything()))),
              ~ case_when(. %like% "1/1|1/0|0/1" ~ "1",
                          . %like% "2/2|2/0|0/2" ~ "2",
                          . %like% "3/3|3/0|0/3" ~ "3",
                          . %like% "4/4|4/0|0/4" ~ "4",
                          . %like% "0/0"          ~ "0",
                          . %like% "./."          ~ "-1",
                          TRUE                    ~ .)) |>
    dplyr::select(-c(3:9)) |>
    arrange(CHROM, POS) |>
    write_tsv(out_tsv)
  invisible(out_tsv)
}

# Convenience: autodetect the number of ## header lines in a plain-text VCF.
count_vcf_header <- function(vcf_path) {
  n <- 0L
  con <- file(vcf_path, "r")
  on.exit(close(con))
  repeat {
    ln <- readLines(con, n = 1L)
    if (length(ln) == 0 || !startsWith(ln, "##")) break
    n <- n + 1L
  }
  n
}

# For cluster-specific: iterate over *_grep_patterns.tsv files.
make_cluster_hmmibd_tsvs <- function(vcf_path, cluster_dir, pattern = "_grep_patterns.tsv$") {
  files <- list.files(cluster_dir, pattern = pattern, full.names = TRUE)
  for (f in files) {
    out <- str_replace(f, "_grep_patterns.tsv$", "_hmmIBD.tsv")
    message(sprintf("[17] %s -> %s", basename(f), basename(out)))
    make_hmmibd_tsv(vcf_path, f, out)
  }
  invisible(files)
}
