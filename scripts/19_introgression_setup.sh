#!/usr/bin/env bash
# Stage B-9 setup: build the hmmIBD.tsv genotype file for the introgression
# analysis. Filters no_MOI VCF to the clusters_combined SNPs (union of the
# per-cluster --maf 0.05 variants — Jacob's original whole-cohort introgression
# input), decompresses, and hands off to 17_hmmibd_genotype.R.
#
# NB (per prompt): Jacob is rewriting the introgression script to be data-
# agnostic, so this pass is a faithful lift — don't over-invest in exact matching.
set -euo pipefail

cd "$(dirname "$0")/.."
CM=data/processed/cluster_maf
NO_MOI=data/gadi/Consensus_SNPs_subset_no_MOI.vcf.gz
OUT=data/processed/introgression
mkdir -p "$OUT"

command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH" >&2; exit 1; }
[[ -s "$CM/clusters_combined.bim" ]] || { echo "ERROR: run 12_cluster_maf.sh first" >&2; exit 1; }

# ---- grep_patterns.tsv from clusters_combined.bim ---------------------------
# BIM has numeric chroms (1..14) after stage-2 rename; VCF has ordered_PKNH_NN_v2.
awk 'BEGIN{OFS="\t"} {printf "ordered_PKNH_%02d_v2\t%s\n", $1, $4}' \
    "$CM/clusters_combined.bim" > "$OUT/grep_patterns.tsv"
echo "grep_patterns lines = $(wc -l < $OUT/grep_patterns.tsv)"

# ---- filtered VCF (samples: all of no_MOI; SNPs: clusters_combined union) ---
if [[ -s "$OUT/Consensus_SNPs_subset_no_MOI.vcf" ]]; then
  echo "--- filtered VCF exists — skipping ---"
else
  echo "--- bcftools view (SNPs only, plain text VCF) ---"
  bcftools view -R "$OUT/grep_patterns.tsv" -Ov "$NO_MOI" \
      > "$OUT/Consensus_SNPs_subset_no_MOI.vcf"
fi
echo "filtered VCF size: $(du -h $OUT/Consensus_SNPs_subset_no_MOI.vcf | cut -f1)"
head -1 "$OUT/Consensus_SNPs_subset_no_MOI.vcf" | head -c 80; echo "..."
head -100 "$OUT/Consensus_SNPs_subset_no_MOI.vcf" | grep -c "^##" > "$OUT/n_header_lines.txt" || true
echo "n_header_lines (##) = $(cat $OUT/n_header_lines.txt)"

# ---- run VCF → hmmIBD.tsv via 17_hmmibd_genotype.R --------------------------
echo "--- Rscript 17_hmmibd_genotype.R (make_hmmibd_tsv) ---"
Rscript -e '
  source("scripts/_setup.R")
  source("scripts/17_hmmibd_genotype.R")
  vcf <- file.path("data/processed/introgression", "Consensus_SNPs_subset_no_MOI.vcf")
  n   <- count_vcf_header(vcf)
  cat(sprintf("[19] autodetected %d header lines\n", n))
  make_hmmibd_tsv(vcf,
                  file.path("data/processed/introgression", "grep_patterns.tsv"),
                  file.path("data/processed/introgression", "hmmIBD.tsv"),
                  skip_header_lines = n)
  cat("[19] hmmIBD.tsv rows:", length(readLines(file.path("data/processed/introgression","hmmIBD.tsv"))), "\n")
'
echo "--- DONE ---"
ls -la "$OUT"
