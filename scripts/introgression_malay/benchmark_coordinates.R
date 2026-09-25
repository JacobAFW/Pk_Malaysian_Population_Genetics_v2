#!/usr/bin/env Rscript
# benchmark_coordinates.R — score this run's filtered windows against the
# published (poster) introgression windows BY GENOMIC COORDINATE.
# --------------------------------------------------------------------------
# WHY BY COORDINATE: the published `WINDOW` ids are row numbers over the sorted
# (CHROM, bin) pairs, so they shift whenever the variant set shifts. Comparing
# ids is what produced the "19% overlap" artifact (agnostic spec, Method status
# block); compared by CHROM + window start the same density method reproduces
# 138/217 = 64% of its own published windows. This script never touches ids.
#
# GRID ALIGNMENT (the one silent-failure mode, checked at runtime):
#   truth:  start = 40000, end = 49999          -> start is ON the 10 kb grid
#   ours:   BIN = 135000, START = 130000, END = 139999
#   BIN is the bin MIDPOINT (START + window/2). Matching truth$start to BIN
#   would compare points 5 kb apart and score ~0. We match truth$start to our
#   START, and assert both are congruent to 0 mod window_size before scoring.
#
# CLI:
#   Rscript benchmark_coordinates.R --config config/malay_cohort.yaml \
#     --out results/<run>/benchmark/coordinate_benchmark.tsv
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(yaml)
})
setDTthreads(3L)

parse_args <- function(av) {
  out <- list(); i <- 1
  while (i <= length(av)) { out[[sub("^--", "", av[i])]] <- av[i + 1]; i <- i + 2 }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
for (r in c("config", "out")) if (is.null(args[[r]])) stop("Missing --", r)
cfg <- yaml::read_yaml(args[["config"]])
RES <- cfg$outputs$results_dir
WIN <- as.integer(cfg$introgression$window_size_bp)
dir.create(dirname(args[["out"]]), recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(sprintf("[benchmark] %s", sprintf(...)))

# NB: sprintf with %d, never paste0. R renders round doubles in scientific
# notation (as.character(100000) == "1e+05"), so a paste0 key built from a
# double silently fails to match the same coordinate read as an integer.
key <- function(chrom, start) sprintf("%d:%d", as.integer(chrom), as.integer(start))

# -- truth -----------------------------------------------------------------
truth <- read_tsv(file.path(RES, "inputs/truth/truth_windows_coords.tsv"),
                  show_col_types = FALSE) %>%
  mutate(CHROM = as.integer(CHROM), start = as.numeric(start))
stopifnot(all(truth$start %% WIN == 0))
say("truth: %d windows (%s)", nrow(truth),
    paste(sprintf("%s %d", names(table(truth$truth_cluster)),
                  as.integer(table(truth$truth_cluster))), collapse = ", "))

# -- eligibility universe: which truth windows this cohort could reach at all
# Same definition the floor derivation uses: a (sample, window) is eligible
# where the sample has more than min_snps non-missing calls. A truth window
# outside this universe cannot be reproduced for reasons that have nothing to
# do with the detection method.
say("recomputing the eligible window universe from the genotype table ...")
gt <- fread(path.expand(cfg$inputs$genotype), header = TRUE, sep = "\t")
setnames(gt, 1:2, c("CHROM", "POS"))
clus <- read_tsv(file.path(RES, "inputs/admix_clusters.tsv"), show_col_types = FALSE)
samp <- intersect(clus$Sample, names(gt))
min_snps <- as.integer(cfg$introgression$min_snps_per_window)
gt[, START := (as.integer(POS) %/% WIN) * WIN]
elig <- gt[, {
  nn <- vapply(.SD, function(v) sum(v != -1L), integer(1))
  .(n_samples_eligible = sum(nn > min_snps))
}, by = .(CHROM = as.integer(CHROM), START), .SDcols = samp][n_samples_eligible > 0]
universe <- key(elig$CHROM, elig$START)
say("eligible universe: %d windows", length(universe))
rm(gt); invisible(gc())

truth <- truth %>% mutate(k = key(CHROM, start), attainable = k %in% universe)
say("truth windows inside the eligible universe: %d of %d",
    sum(truth$attainable), nrow(truth))

# -- our called windows, per arm and per cluster ---------------------------
read_called <- function(arm) {
  f <- file.path(RES, "aggregate", arm, "introgressed_windows_filtered.tsv")
  d <- fread(f, sep = "\t", select = c("Cluster", "CHROM", "START"))
  unique(d[, .(Cluster, k = key(CHROM, START))])
}

score <- function(called_k, truth_sub, label_arm, label_dir, scope) {
  tk <- truth_sub$k
  shared <- length(intersect(called_k, tk))
  union_n <- length(union(called_k, tk))
  tibble::tibble(
    arm = label_arm, direction = label_dir, scope = scope,
    truth_windows = length(tk), called_windows = length(called_k),
    reproduced = shared,
    pct_reproduced = ifelse(length(tk) > 0, 100 * shared / length(tk), NA_real_),
    jaccard = ifelse(union_n > 0, shared / union_n, NA_real_)
  )
}

rows <- list()
for (arm in c("full", "declonal")) {
  called <- read_called(arm)
  all_k <- unique(called$k)
  for (scope in c("all_truth", "attainable_only")) {
    tsub_all <- if (scope == "all_truth") truth else truth %>% filter(attainable)
    # combined: our whole filtered window set vs the whole truth set
    rows[[length(rows) + 1]] <- score(all_k, tsub_all, arm, "combined", scope)
    # per truth direction, scored against OUR matching cluster's windows
    for (cl in c("Mf", "Mn")) {
      ck <- called[Cluster == cl, k]
      rows[[length(rows) + 1]] <- score(ck, tsub_all %>% filter(truth_cluster == cl),
                                        arm, cl, scope)
    }
    # per truth direction, scored against our FULL window set (a truth window
    # may be recovered under a different cluster label than the poster gave it)
    for (cl in c("Mf", "Mn")) {
      rows[[length(rows) + 1]] <- score(all_k, tsub_all %>% filter(truth_cluster == cl),
                                        arm, paste0(cl, "_vs_anycluster"), scope)
    }
  }
}
# -- supplementary: the SAME scoring at the published chain's own per-cluster
# floor (n > 5, i.e. N = 6), from the phase-3 sweep work dirs. This is not a
# re-tune of the headline — the derived floor stands — it exists to separate
# "our floor is stricter" from "our windows are elsewhere" when reading the
# gap to the published result. Clearly labelled as supplementary.
for (arm in c("full", "declonal")) {
  p <- file.path(RES, "floor", paste0("_sweep_", arm), "N06",
                 "introgressed_windows_filtered.tsv")
  if (!file.exists(p)) { say("supplementary N=6 table missing for %s — skipped", arm); next }
  d <- fread(p, sep = "\t", select = c("Cluster", "CHROM", "START"))
  d <- unique(d[, .(Cluster, k = key(CHROM, START))])
  rows[[length(rows) + 1]] <- score(unique(d$k), truth, arm, "combined",
                                    "published_floor_6_supplementary")
  for (cl in c("Mf", "Mn")) {
    rows[[length(rows) + 1]] <- score(d[Cluster == cl, k],
                                      truth %>% filter(truth_cluster == cl),
                                      arm, cl, "published_floor_6_supplementary")
  }
}

out <- bind_rows(rows)
write_tsv(out, args[["out"]])

# -- context table: what the comparison is against ------------------------
fam <- cfg$truth$mn_windows  # same fixture dir as the copied truth
fam <- file.path(dirname(dirname(path.expand(fam))), "inputs", "cleaned.fam")
ctx <- tibble::tibble(item = character(), value = character())
add <- function(i, v) ctx <<- bind_rows(ctx, tibble::tibble(item = i, value = as.character(v)))
add("our_cluster_basis_n", nrow(clus))
if (file.exists(fam)) {
  f558 <- unique(read.table(fam, header = FALSE, stringsAsFactors = FALSE)[[1]])
  add("published_basis_n", length(f558))
  add("basis_intersection", length(intersect(f558, clus$Sample)))
  add("in_published_not_ours", length(setdiff(f558, clus$Sample)))
  add("in_ours_not_published", length(setdiff(clus$Sample, f558)))
}
add("truth_windows_total", nrow(truth))
add("truth_windows_attainable_here", sum(truth$attainable))
add("eligible_universe_windows", length(universe))
add("published_558basis_ceiling", "138/217 = 64% (agnostic spec, regen_density)")
write_tsv(ctx, file.path(dirname(args[["out"]]), "benchmark_context.tsv"))

print(as.data.frame(out %>% filter(scope != "attainable_only")), row.names = FALSE)
print(as.data.frame(ctx), row.names = FALSE)
say("wrote %s (+ benchmark_context.tsv)", args[["out"]])
