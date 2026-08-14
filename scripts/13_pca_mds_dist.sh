#!/usr/bin/env bash
# Stage B-5: PCA + MDS + IBS distance on the poster cohort (Analyses.Rmd L770-773).
# Uses cleaned.{bed,bim,fam} (regen 559; F-anchored labels applied at plot time).
set -euo pipefail

cd "$(dirname "$0")/.."
PLINK=data/processed/plink
OUT=data/processed/ordination
mkdir -p "$OUT"

command -v plink >/dev/null || { echo "ERROR: plink not on PATH" >&2; exit 1; }

echo "--- PLINK --cluster --mds-plot 10 ---"
plink --bfile "$PLINK/cleaned" --cluster --mds-plot 10 --allow-no-sex --out "$OUT/Pk"

echo "--- PLINK --pca ---"
plink --bfile "$PLINK/cleaned" --pca --allow-no-sex --out "$OUT/Pk"

echo "--- PLINK --distance square ---"
plink --bfile "$PLINK/cleaned" --distance square --allow-no-sex --out "$OUT/Pk"

echo "--- DONE ---"
ls -la "$OUT"
