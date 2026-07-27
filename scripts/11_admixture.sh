#!/usr/bin/env bash
# Stage B-3: ADMIXTURE K=1..6 with CV, on cleaned.bed (poster cohort = 558 samples).
# Inputs : data/processed/plink/cleaned.{bed,bim,fam}
# Outputs: data/processed/admixture/cleaned.{K}.{P,Q}   for K in 1..6
#          data/processed/admixture/log{K}.out
#          data/processed/admixture/Pk.cv.error          (K, CV) two-col
# Run    : source vvg-box-pixi/bin/activate; bash scripts/11_admixture.sh
set -euo pipefail

cd "$(dirname "$0")/.."
PLINK=data/processed/plink
OUT=data/processed/admixture
mkdir -p "$OUT"

[[ -s "$PLINK/cleaned.bed" ]] || { echo "ERROR: $PLINK/cleaned.bed missing; run 10_normalise_plink.sh first" >&2; exit 1; }
command -v admixture >/dev/null || { echo "ERROR: admixture not on PATH" >&2; exit 1; }
admixture --version 2>&1 | head -1 || true

# ADMIXTURE writes outputs (cleaned.K.{P,Q}) to its CWD; run from $OUT, point at the bed.
# From data/processed/admixture/, the bed is at ../plink/cleaned.bed.
BED_REL="../plink/cleaned.bed"

# Ensure chrom column 1 = "0"-equivalent (ADMIXTURE rejects non-human chrom names).
# Our bim has integer chroms 1..14 already (from chrom rename in stage B-2); that's fine
# for ADMIXTURE (it accepts 0..22 numeric).
head -1 "$PLINK/cleaned.bim" | awk '{print "first bim row chrom = "$1}'

pushd "$OUT" >/dev/null
for K in 1 2 3 4 5 6; do
  log=log${K}.out
  if [[ -s cleaned.${K}.Q && -s "$log" ]]; then
    echo "--- K=$K: outputs exist (skipping); rm cleaned.${K}.* to redo ---"
  else
    echo "--- ADMIXTURE K=$K (cv) → $log ---"
    admixture --cv "$BED_REL" $K > "$log"
  fi
done

# Collate CV errors
awk '/CV/ {print $3, $4}' log*.out | sed 's/[()K=:]//g' | awk '{print $1"\t"$2}' \
  > Pk.cv.error
echo "--- CV table (lower = better) ---"
sort -n -k1,1 Pk.cv.error
popd >/dev/null

echo "--- DONE ---"
ls -la "$OUT"
