#!/usr/bin/env bash
# Stage B-7: Fst Mn vs Mf on clusters_combined.
# Uses truth-anchored cluster labels (F decision) from data/processed/cluster_maf/
# cluster_labels.tsv. Peninsular is excluded from the two-population Fst.
set -euo pipefail

cd "$(dirname "$0")/.."
CM=data/processed/cluster_maf
OUT=data/processed/fst
mkdir -p "$OUT"

command -v plink >/dev/null || { echo "ERROR: plink not on PATH" >&2; exit 1; }
[[ -s "$CM/clusters_combined.bed" ]] || { echo "ERROR: run 12_cluster_maf.sh first" >&2; exit 1; }
[[ -s "$CM/cluster_labels.tsv" ]] || { echo "ERROR: cluster_labels.tsv missing" >&2; exit 1; }

# ---- sample_clusters.txt: FID IID Cluster (all 501 samples) -----------------
awk -F'\t' 'NR>1 {print $1"\t"$1"\t"$2}' "$CM/cluster_labels.tsv" > "$OUT/sample_clusters.txt"
echo "sample_clusters lines = $(wc -l < $OUT/sample_clusters.txt)"

# ---- exclude Peninsular for the two-pop Fst ---------------------------------
awk -F'\t' 'NR>1 && $2=="Peninsular" {print $1"\t"$1}' "$CM/cluster_labels.tsv" > "$OUT/exclude_pen.txt"
echo "exclude_pen lines     = $(wc -l < $OUT/exclude_pen.txt)"

echo "--- plink --fst (Mn vs Mf) ---"
plink --bfile "$CM/clusters_combined" \
      --within "$OUT/sample_clusters.txt" \
      --remove "$OUT/exclude_pen.txt" \
      --fst \
      --allow-no-sex \
      --out "$OUT/Pk"

echo "--- DONE ---"
head -3 "$OUT/Pk.fst"
echo "..."
echo "Fst rows (per-SNP) = $(wc -l < $OUT/Pk.fst)"
echo "plink summary:"
grep -E "Fst estimate|weighted" "$OUT/Pk.log"
