#!/usr/bin/env bash
# Phase 5 / FEEMS step 1 driver: build the FEEMS genotype input.
#
#   2x_feems_prep.py   -> keep-list (cohort INTERSECT EEMS keep-list INTERSECT region), coords, outer ring
#   plink --keep       -> restrict to that sample set
#   plink --geno/--maf -> drop high-missingness SNPs, then MONOMORPHIC SNPs  <-- the feems #30 fix
#
# THE #30 FIX: feems scales each SNP by 1/sqrt(mu*(1-mu)). A monomorphic site has
# mu == 0, so the scaling divides by zero -> inf/nan in the sample covariance ->
# "did not converge" AssertionError inside L-BFGS. Monomorphic sites MUST be gone
# before SpatialGraph is constructed. `--maf $MAF` (default 1e-6) does that here;
# 2x_feems_fit.py re-asserts zero-variance-free on the loaded matrix.
#
# Usage:
#   bash scripts/2x_feems_prep.sh                       # primary  (cleaned_501 source)
#   SOURCE=clusters_combined bash scripts/2x_feems_prep.sh   # sensitivity run
#
# Env overrides: SOURCE, TAG, MAF, GENO, BUFFER
set -euo pipefail
cd "$(dirname "$0")/.."

PLINK=${PLINK:-vvg-box-pixi/opt/pixi/global/bin/plink}
# micromamba is a shell function in interactive shells; use the binary directly here.
MM=${MM:-${MAMBA_EXE:-$HOME/mamba/micromamba}}
ENV=${ENV:-feems_e}

SOURCE=${SOURCE:-cleaned_501}          # cleaned_501 | clusters_combined
TAG=${TAG:-$SOURCE}
MAF=${MAF:-1e-6}                       # >0 == drop monomorphic sites (the #30 fix)
GENO=${GENO:-0.1}                      # max per-SNP missingness
BUFFER=${BUFFER:-0.25}                 # degrees of buffer on the region polygon
CLUSTER=${CLUSTER:-}                   # optional: restrict to one ancestry cluster (Mf/Mn/Peninsular)
CLUSTER_LABELS=${CLUSTER_LABELS:-data/processed/cluster_maf/cluster_labels.tsv}

BFILE=data/processed/cluster_maf/${SOURCE}
KEEPLIST=data/raw/eems/samples_to_keep.txt
COORDS=data/raw/eems/gis_ordered_for_python.coords
OUTER=data/raw/MSC/spatial/malaysia.shp
OUT=data/processed/feems

[[ -x "$PLINK" ]] || { echo "ERROR: plink not found/executable at $PLINK" >&2; exit 1; }
for f in "$BFILE.bed" "$BFILE.bim" "$BFILE.fam" "$KEEPLIST" "$COORDS" "$OUTER"; do
  [[ -s "$f" ]] || { echo "ERROR: missing input $f" >&2; exit 1; }
done
mkdir -p "$OUT"

echo "=== [1/3] assemble keep-list + coords + outer ring (tag=$TAG, source=$SOURCE) ==="
CLUSTER_ARGS=()
[[ -n "$CLUSTER" ]] && CLUSTER_ARGS=(--cluster "$CLUSTER" --cluster-labels "$CLUSTER_LABELS")
$MM run -n "$ENV" python scripts/2x_feems_prep.py \
    --bfile "$BFILE" --keeplist "$KEEPLIST" --coords "$COORDS" \
    --outer-shp "$OUTER" --buffer "$BUFFER" --out-dir "$OUT" --tag "$TAG" "${CLUSTER_ARGS[@]}"

echo
echo "=== [2/3] plink subset -> --geno $GENO -> --maf $MAF (drops monomorphic sites) ==="
$PLINK --bfile "$BFILE" \
       --keep "$OUT/${TAG}_keep.txt" \
       --geno "$GENO" \
       --maf "$MAF" \
       --make-bed --allow-no-sex \
       --out "$OUT/${TAG}_feems" >/dev/null
grep -E "variants|people|genotyping rate|removed due to" "$OUT/${TAG}_feems.log" | sed 's/^/  /'

echo
echo "=== [3/3] summary ==="
printf '  source bfile      : %s (%s variants, %s samples)\n' \
  "$BFILE" "$(wc -l < "$BFILE.bim" | tr -d ' ')" "$(wc -l < "$BFILE.fam" | tr -d ' ')"
printf '  feems bfile       : %s_feems (%s variants, %s samples)\n' \
  "$OUT/$TAG" "$(wc -l < "$OUT/${TAG}_feems.bim" | tr -d ' ')" "$(wc -l < "$OUT/${TAG}_feems.fam" | tr -d ' ')"
printf '  filters           : --geno %s --maf %s   (--maf > 0 == monomorphic sites dropped)\n' "$GENO" "$MAF"
echo "  next: micromamba run -n $ENV python scripts/2x_feems_fit.py --tag $TAG"
