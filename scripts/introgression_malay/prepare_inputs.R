#!/usr/bin/env Rscript
# prepare_inputs.R — build the derived input tables the lifted introgression
# module expects, from this cohort's own files.
# --------------------------------------------------------------------------
# WHAT: turns the project's raw/processed inputs into the five tables
#       scripts/introgression_module/*.R read, plus the de-clonalization
#       drop-list and the coordinate-truth copy used by the benchmark.
# WHY:  the module is cohort-agnostic and takes generic table shapes; every
#       cohort needs a thin prep step. Nothing here reimplements module logic.
#
# Writes ONLY under --out-dir (a gitignored results/<run>/inputs path).
# Reads data/ read-only.
#
# CLI:
#   Rscript prepare_inputs.R --config config/malay_cohort.yaml \
#                            --out-dir results/<run>/inputs
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(yaml)
})

setDTthreads(3L)  # protocol: at most 3 threads per tool

parse_args <- function(av) {
  out <- list(); i <- 1
  while (i <= length(av)) {
    key <- av[i]
    if (!startsWith(key, "--")) stop("Bad arg: ", key)
    out[[sub("^--", "", key)]] <- av[i + 1]; i <- i + 2
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
for (req in c("config", "out-dir")) if (is.null(args[[req]])) stop("Missing --", req)
OUT <- args[["out-dir"]]
dir.create(file.path(OUT, "truth"), recursive = TRUE, showWarnings = FALSE)

say <- function(...) message(sprintf("[prepare_inputs] %s", sprintf(...)))

# Every cohort-specific path comes from the config (config/malay_cohort.yaml).
# Nothing about this cohort — reference, annotation, genotype, labels, truth —
# is written into the code.
cfg <- yaml::read_yaml(args[["config"]])
P <- c(
  cfg$inputs,
  list(truth_filt = cfg$truth$filtered,
       truth_mf   = cfg$truth$mf_windows,
       truth_mn_src = cfg$truth$mn_windows)
)
P <- lapply(P, path.expand)
for (nm in names(P)) if (!file.exists(P[[nm]])) stop("Input not found (", nm, "): ", P[[nm]])
say("config: %s (run_id %s)", args[["config"]], cfg$run_id)

# --------------------------------------------------------------------------
# 1. Clusters + metadata
# --------------------------------------------------------------------------
labels <- read_tsv(P$labels, show_col_types = FALSE)
stopifnot(all(c("Sample", "Cluster") %in% names(labels)))
say("labels: %d samples — %s", nrow(labels),
    paste(sprintf("%s n=%d", names(table(labels$Cluster)),
                  as.integer(table(labels$Cluster))), collapse = ", "))

gt_header <- names(fread(P$genotype, nrows = 0, sep = "\t"))
gt_samples <- gt_header[-(1:2)]
present <- labels$Sample %in% gt_samples
if (!all(present)) stop(sum(!present), " label samples absent from the genotype table")
say("genotype table: %d sample columns; all %d labelled samples present",
    length(gt_samples), nrow(labels))

clusters <- labels %>% dplyr::select(Sample, Cluster) %>% arrange(Sample)
write_tsv(clusters, file.path(OUT, "admix_clusters.tsv"))
# Minimal metadata: sample_id only. No `geography` role, so the aggregate
# step's per-geography summary skips itself with a logged note (by design).
write_tsv(tibble::tibble(sample_id = clusters$Sample), file.path(OUT, "samples.tsv"))

# --------------------------------------------------------------------------
# 2. Contig map (headerless CONTIG -> integer CHROM) + .fai reference
# --------------------------------------------------------------------------
# The GFF and the reference FASTA use the same contig names here
# (ordered_PKNH_NN_v2), so the aggregate step's mapping resolves by exact name
# and never needs the length fallbacks. CHROM codes are fai order over the
# nuclear contigs, 1-based — which is what hmmIBD.tsv's "01".."14" encode.
fai <- fread(P$fai, header = FALSE, sep = "\t", select = 1:2,
             col.names = c("CONTIG", "LEN"))
nuclear <- fai[!grepl("MIT|API", CONTIG)]
cmap <- nuclear[, .(CONTIG, CHROM = seq_len(.N))]
fwrite(cmap, file.path(OUT, "contig_map.tsv"), sep = "\t", col.names = FALSE)
say("contig_map: %d nuclear contigs (%d non-nuclear dropped)",
    nrow(cmap), nrow(fai) - nrow(nuclear))

gt_chrom <- unique(fread(P$genotype, sep = "\t", select = 1)[[1]])
stopifnot(setequal(as.integer(gt_chrom), cmap$CHROM))
say("CHROM codes in the genotype table match contig_map exactly (%d contigs)", nrow(cmap))

# --------------------------------------------------------------------------
# 3. Per-sample missingness (the `.smiss` substitute — see the phase report)
# --------------------------------------------------------------------------
# The representative rule in the module's docs/clonality.md is "lowest
# per-sample missingness (PLINK .smiss), then alphabetical". No .smiss/.imiss
# exists anywhere in this project, so F_MISS is computed directly from the
# genotype table (missing = -1) — the same quantity PLINK would report, from
# the same calls, deterministically. The documented degradation is
# alphabetical-only; both choices are recorded so the difference is measurable.
say("computing per-sample missingness from the genotype table ...")
t0 <- proc.time()[["elapsed"]]
gt <- fread(P$genotype, header = TRUE, sep = "\t")
setnames(gt, 1:2, c("CHROM", "POS"))
n_sites <- nrow(gt)
miss <- vapply(clusters$Sample, function(s) sum(gt[[s]] == -1L), numeric(1))
smiss <- tibble::tibble(Sample = clusters$Sample,
                        N_MISS = as.integer(miss),
                        F_MISS = miss / n_sites) %>% arrange(Sample)
write_tsv(smiss, file.path(OUT, "sample_missingness.tsv"))
say("missingness over %d sites: median F_MISS %.4f, range %.4f-%.4f (%.1fs)",
    n_sites, median(smiss$F_MISS), min(smiss$F_MISS), max(smiss$F_MISS),
    proc.time()[["elapsed"]] - t0)
rm(gt); invisible(gc())

# --------------------------------------------------------------------------
# 4. De-clonalization: connected components -> one representative per group
# --------------------------------------------------------------------------
read_pairs <- function(path) {
  d <- fread(path, header = TRUE, sep = "\t", select = 1:2,
             col.names = c("s1", "s2"), colClasses = "character")
  unique(d[!is.na(s1) & !is.na(s2) & s1 != s2])
}

# The two clonal pair tables MISSION.md §3 names are in the field-collection id
# namespace, not the sequencing-id namespace the genotype table and the labels
# use. Checked rather than assumed, and reported either way.
named_pairs <- unique(rbindlist(list(read_pairs(P$clonal_all), read_pairs(P$clonal_mn))))
named_ids <- sort(unique(c(named_pairs$s1, named_pairs$s2)))
n_named_usable <- length(intersect(named_ids, clusters$Sample))
say("named clonal tables: %d pairs over %d samples; %d of those ids are in the clustered set",
    nrow(named_pairs), length(named_ids), n_named_usable)

# De-clonalization source: this project's full pairwise hmmIBD output, which is
# in the genotype id namespace. Threshold comes from the config, which cites
# this project's own rule (scripts/14b_clonal_networks.R).
thr <- as.numeric(cfg$clonality$threshold)
ibd_col <- cfg$clonality$ibd_column
ibd_path <- path.expand(cfg$clonality$pairwise_ibd)
if (!file.exists(ibd_path)) stop("Pairwise IBD table not found: ", ibd_path)
ibd <- fread(ibd_path, header = TRUE, sep = "\t",
             select = c(cfg$clonality$id_columns, ibd_col))
setnames(ibd, c("s1", "s2", "ibd"))
say("pairwise IBD: %d pairs over %d samples (%s)", nrow(ibd),
    length(unique(c(ibd$s1, ibd$s2))), basename(ibd_path))
pairs <- unique(ibd[ibd >= thr & s1 != s2,
                    .(s1 = as.character(s1), s2 = as.character(s2))])
n_pairs_raw <- nrow(pairs)
# Only samples that are actually in this analysis can be collapsed.
pairs <- pairs[s1 %in% clusters$Sample & s2 %in% clusters$Sample]
say("clonal pairs at %s >= %.3f: %d total, %d with both members in the %d-sample clustered set",
    ibd_col, thr, n_pairs_raw, nrow(pairs), nrow(clusters))
# Threshold separation: how many pairs sit just below the cut. A crowded
# shoulder would mean the cut is arbitrary; a gap means it is not.
near <- nrow(ibd[ibd >= thr - 0.08 & ibd < thr])
say("pairs in [%.2f, %.2f): %d — threshold separation check", thr - 0.08, thr, near)
if (nrow(pairs) == 0) stop("No clonal pairs within the clustered set — cannot build the de-clonalized arm")

# Union-find over the retained pairs.
parent <- new.env(parent = emptyenv())
find <- function(x) {
  if (is.null(parent[[x]])) { parent[[x]] <- x; return(x) }
  while (parent[[x]] != x) { parent[[x]] <- parent[[parent[[x]]]]; x <- parent[[x]] }
  x
}
unite <- function(a, b) { ra <- find(a); rb <- find(b); if (ra != rb) parent[[ra]] <- rb }
for (i in seq_len(nrow(pairs))) unite(pairs$s1[i], pairs$s2[i])

members <- sort(unique(c(pairs$s1, pairs$s2)))
audit <- tibble::tibble(Sample = members,
                        root = vapply(members, find, character(1))) %>%
  left_join(clusters, by = "Sample") %>%
  left_join(smiss %>% dplyr::select(Sample, F_MISS), by = "Sample")

# Group id = alphabetically-first member, so it is stable across re-runs
# regardless of union-find's internal root choice.
audit <- audit %>% group_by(root) %>% mutate(group_id = min(Sample)) %>% ungroup()

# Straddling groups: representative comes from the group's majority cluster;
# a tie between clusters goes to the alphabetically-first cluster name.
group_cluster <- audit %>%
  count(group_id, Cluster, name = "n_in_cluster") %>%
  group_by(group_id) %>%
  arrange(desc(n_in_cluster), Cluster, .by_group = TRUE) %>%
  summarise(majority_cluster = first(Cluster),
            n_clusters = dplyr::n(), .groups = "drop")
n_straddle <- sum(group_cluster$n_clusters > 1)
if (n_straddle > 0) say("WARNING: %d clonal group(s) straddle more than one cluster", n_straddle)

# Representative: lowest F_MISS within the majority cluster, alphabetical tie-break.
reps <- audit %>%
  left_join(group_cluster, by = "group_id") %>%
  group_by(group_id) %>%
  arrange(Cluster != majority_cluster, F_MISS, Sample, .by_group = TRUE) %>%
  mutate(is_representative = dplyr::row_number() == 1L,
         # what alphabetical-only (the documented no-.smiss fallback) would pick
         rep_alphabetical = Sample == min(Sample[Cluster == majority_cluster])) %>%
  ungroup()

audit_out <- reps %>%
  dplyr::select(group_id, Sample, Cluster, majority_cluster, n_clusters,
                F_MISS, is_representative, rep_alphabetical) %>%
  arrange(group_id, desc(is_representative), Sample)
write_tsv(audit_out, file.path(OUT, "declonalization_audit.tsv"))

drop <- audit_out %>% filter(!is_representative) %>% pull(Sample) %>% sort()
writeLines(drop, file.path(OUT, "exclude_unique.txt"))
writeLines(character(0), file.path(OUT, "exclude_full.txt"))
writeLines(sort(clusters$Sample), file.path(OUT, "all_genotypes.txt"))
writeLines(sort(setdiff(clusters$Sample, drop)), file.path(OUT, "unique_genotypes.txt"))

n_groups <- dplyr::n_distinct(audit_out$group_id)
n_disagree <- audit_out %>% filter(is_representative != rep_alphabetical) %>%
  dplyr::distinct(group_id) %>% nrow()
say("clonal groups: %d covering %d samples; %d dropped, %d kept as representatives",
    n_groups, nrow(audit_out), length(drop), n_groups)
say("representative choice differs from alphabetical-only in %d of %d groups",
    n_disagree, n_groups)

collapse <- clusters %>%
  mutate(kept = !(Sample %in% drop)) %>%
  group_by(Cluster) %>%
  summarise(n_full = dplyr::n(), n_declonal = sum(kept),
            n_dropped = sum(!kept), .groups = "drop") %>%
  arrange(desc(n_full))
write_tsv(collapse, file.path(OUT, "declonalization_counts.tsv"))
print(as.data.frame(collapse))

# --------------------------------------------------------------------------
# 5. Coordinate truth for the benchmark (CHROM + window start, never WINDOW id)
# --------------------------------------------------------------------------
file.copy(P$truth_mf, file.path(OUT, "truth", "mf_windows.tsv"), overwrite = TRUE)
file.copy(P$truth_mn_src, file.path(OUT, "truth", "mn_windows.tsv"), overwrite = TRUE)
file.copy(P$truth_filt, file.path(OUT, "truth", "introgressed_windows_filtered.tsv"),
          overwrite = TRUE)

tmf <- read_tsv(file.path(OUT, "truth", "mf_windows.tsv"), show_col_types = FALSE)
tmn <- read_tsv(file.path(OUT, "truth", "mn_windows.tsv"), show_col_types = FALSE)
tfl <- read_tsv(file.path(OUT, "truth", "introgressed_windows_filtered.tsv"),
                show_col_types = FALSE)
truth_coords <- bind_rows(tmf %>% mutate(truth_cluster = "Mf"),
                          tmn %>% mutate(truth_cluster = "Mn")) %>%
  mutate(CHROM = as.integer(CHROM), start = as.numeric(start)) %>%
  arrange(CHROM, start)
write_tsv(truth_coords, file.path(OUT, "truth", "truth_windows_coords.tsv"))

ids_filtered <- sort(unique(tfl$WINDOW))
ids_coords   <- sort(unique(c(tmf$WINDOW, tmn$WINDOW)))
say("truth: %d rows / %d unique WINDOW ids in the filtered table; %d Mf + %d Mn coordinate windows",
    nrow(tfl), length(ids_filtered), nrow(tmf), nrow(tmn))
say("coordinate coverage of the filtered WINDOW ids: %d of %d; overlap Mf n Mn = %d",
    length(intersect(ids_filtered, ids_coords)), length(ids_filtered),
    length(intersect(tmf$WINDOW, tmn$WINDOW)))
if (!setequal(ids_filtered, ids_coords))
  say("WARNING: coordinate truth does not cover the filtered WINDOW id set exactly")

# --------------------------------------------------------------------------
# 6. Provenance
# --------------------------------------------------------------------------
md5 <- function(p) {
  o <- suppressWarnings(system2("md5", c("-q", shQuote(p)), stdout = TRUE, stderr = FALSE))
  if (length(o) == 0) NA_character_ else o[1]
}
P$pairwise_ibd <- ibd_path
prov <- tibble::tibble(
  role = c(names(P), "OUT:admix_clusters", "OUT:samples", "OUT:contig_map",
           "OUT:sample_missingness", "OUT:exclude_unique", "OUT:declonalization_audit",
           "OUT:truth_windows_coords"),
  path = c(unlist(P, use.names = FALSE),
           file.path(OUT, c("admix_clusters.tsv", "samples.tsv", "contig_map.tsv",
                            "sample_missingness.tsv", "exclude_unique.txt",
                            "declonalization_audit.tsv", "truth/truth_windows_coords.tsv")))
) %>%
  mutate(bytes = file.size(path), md5 = vapply(path, md5, character(1)),
         mtime = format(file.mtime(path), "%Y-%m-%d %H:%M:%S"))
write_tsv(prov, file.path(OUT, "input_provenance.tsv"))
say("wrote %d provenance rows -> %s", nrow(prov), file.path(OUT, "input_provenance.tsv"))
say("done")
