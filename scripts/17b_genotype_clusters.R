# Pipeline-completeness (not executed for validation this pass):
# Cluster-specific VCF → per-cluster hmmIBD.tsv genotype files.
# Ports scripts/gadi/genotype_clusters.R, same logic + local paths.
#
# Wraps 17_hmmibd_genotype.R's make_hmmibd_tsv() over the per-cluster
# *_grep_patterns.tsv files typically produced by hmmIBD_cluster_spec.pbs.
#
# Usage (bash):
#   Rscript scripts/17b_genotype_clusters.R \
#     <vcf_path> <cluster_dir> [<pattern>]
#
# Example (against a decompressed VCF; produce Mf/Mn/Pen_hmmIBD.tsv):
#   Rscript scripts/17b_genotype_clusters.R \
#     data/processed/hmmibd/Consensus_SNPs_subset_no_MOI.vcf \
#     data/processed/cluster_maf
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts/_setup.R"))
source(here::here("scripts/17_hmmibd_genotype.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: 17b_genotype_clusters.R <vcf_path> <cluster_dir> [<pattern>]")
}
vcf_path   <- args[1]
cluster_dir <- args[2]
pattern    <- if (length(args) >= 3) args[3] else "_grep_patterns.tsv$"

message(sprintf("[17b] VCF: %s", vcf_path))
message(sprintf("[17b] cluster dir: %s  pattern: %s", cluster_dir, pattern))
make_cluster_hmmibd_tsvs(vcf_path, cluster_dir, pattern)
message("[17b] done")
