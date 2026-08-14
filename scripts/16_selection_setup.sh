#!/usr/bin/env bash
# Stage B-8 setup: per-cluster VCFs split by chromosome for rehh iHS.
# Lifted from scripts/gadi/selection_cluster_spec.pbs.
#
# For each cluster (Mf, Mn, Peninsular):
#   1. Build snp_filter.txt = SNPs kept by that cluster's post-MAF .bim
#   2. Build exclude.txt = samples NOT in that cluster (from cluster_labels.tsv)
#   3. bcftools view no_MOI VCF → filter to snps + samples → +missing2ref
#   4. Split per-chromosome (14 nuclear contigs).
set -euo pipefail

cd "$(dirname "$0")/.."
CM=data/processed/cluster_maf
NO_MOI=data/gadi/Consensus_SNPs_subset_no_MOI.vcf.gz
OUT=data/processed/selection
mkdir -p "$OUT"

command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH" >&2; exit 1; }
[[ -s "$CM/cluster_labels.tsv" ]] || { echo "ERROR: run 12_cluster_maf.sh first" >&2; exit 1; }

# List of 14 nuclear contigs (matches Analyses.Rmd + rehh.R)
CONTIGS=(ordered_PKNH_01_v2 ordered_PKNH_02_v2 ordered_PKNH_03_v2 ordered_PKNH_04_v2 \
         ordered_PKNH_05_v2 ordered_PKNH_06_v2 ordered_PKNH_07_v2 ordered_PKNH_08_v2 \
         ordered_PKNH_09_v2 ordered_PKNH_10_v2 ordered_PKNH_11_v2 ordered_PKNH_12_v2 \
         ordered_PKNH_13_v2 ordered_PKNH_14_v2)

for CLUSTER in Mf Mn Pen; do
  CDIR=$OUT/$CLUSTER
  mkdir -p "$CDIR"

  # canonical name in labels.tsv is "Peninsular" (Pen is script alias)
  LABEL=$([[ "$CLUSTER" == "Pen" ]] && echo Peninsular || echo "$CLUSTER")

  echo "--- setup $CLUSTER (label=$LABEL) ---"

  # snp_filter.txt from post-MAF .bim (chrom \t pos). Our .bim uses numeric chroms
  # (1..14) after the stage-2 rename; we need to add "ordered_PKNH_NN_v2" prefix
  # to filter the raw VCF (which still has original chrom names).
  awk 'BEGIN{OFS="\t"} {printf "ordered_PKNH_%02d_v2\t%s\n", $1, $4}' \
      "$CM/$CLUSTER.bim" > "$CDIR/snp_filter.txt"
  echo "  snp_filter lines = $(wc -l < $CDIR/snp_filter.txt)"

  # exclude.txt = FID IID of samples NOT in this cluster (from cluster_labels.tsv)
  awk -F'\t' -v c="$LABEL" 'NR>1 && $2!=c {print $1}' \
      "$CM/cluster_labels.tsv" > "$CDIR/exclude.txt"
  echo "  exclude lines    = $(wc -l < $CDIR/exclude.txt) (samples dropped from no_MOI)"

  # Filter no_MOI VCF: keep only cluster's SNPs + samples, fill missing → ref
  if [[ ! -s "$CDIR/selection.vcf.gz" ]]; then
    echo "  bcftools view + missing2ref → $CDIR/selection.vcf.gz"
    bcftools view --force-samples -R "$CDIR/snp_filter.txt" -S "^$CDIR/exclude.txt" \
             "$NO_MOI" \
      | bcftools +missing2ref -Oz -o "$CDIR/selection.vcf.gz"
    bcftools index -f "$CDIR/selection.vcf.gz"
  else
    echo "  $CDIR/selection.vcf.gz exists — skipping"
  fi

  # per-chrom split
  for i in "${CONTIGS[@]}"; do
    OUT_VCF="$CDIR/${i}.vcf.gz"
    if [[ ! -s "$OUT_VCF" ]]; then
      bcftools view -r "$i" -Oz -o "$OUT_VCF" "$CDIR/selection.vcf.gz"
    fi
  done
  echo "  per-chrom VCFs = $(ls $CDIR/ordered_PKNH_*.vcf.gz | wc -l)"
done

echo "--- DONE ---"
ls -la "$OUT"
