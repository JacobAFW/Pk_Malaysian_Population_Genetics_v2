#!/usr/bin/env bash
# Stage B-4: cluster-specific MAF → clusters_combined_fixed.ped
# Lifted from scripts/gadi/{cluster_specific_maf.Rmd, Analyses.Rmd L1223-1267}.
#
# Per Jacob's F decision (2026-07-27): use truth data/raw/tess3r/cleaned.3.Q for
# cluster labels (anchor downstream to poster-faithful membership). Since our
# regen cohort has 559 samples but truth has 558 (only 501 ∩), we restrict to
# the 501 ∩ samples via --keep so every downstream sample has a truth label.
#
# Inputs : data/processed/plink/cleaned.{bed,bim,fam}     (regen 559)
#          data/raw/tess3r/cleaned.{fam,3.Q}              (truth 558 → cluster labels)
#          data/raw/tess3r/complete_gis_na_omit.fam        (GIS-non-missing 532)
# Outputs: data/processed/cluster_maf/
#            - Exclude_for_{Mf,Mn,Pen}.tsv (from truth cluster labels for 501 ∩)
#            - Mf.{bed,bim,fam}, Mn.*, Pen.*
#            - clusters_combined.{bed,bim,fam,ped,map}
#            - clusters_combined_fixed.ped (post sed */0)
# Run    : source vvg-box-pixi/bin/activate; bash scripts/12_cluster_maf.sh
set -euo pipefail

cd "$(dirname "$0")/.."
PLINK=data/processed/plink
OUT=data/processed/cluster_maf
mkdir -p "$OUT"

command -v plink >/dev/null || { echo "ERROR: plink not on PATH; source vvg-box-pixi/bin/activate" >&2; exit 1; }

# ---- (1) derive truth-anchored cluster labels + Exclude files ---------------
echo "--- derive cluster labels from truth cleaned.3.Q (F-anchored) ---"
Rscript scripts/12b_derive_cluster_labels.R

# The R script writes:
#   $OUT/keep_501.txt         — FID IID (2 col) for --keep, the 501 ∩ samples
#   $OUT/cluster_labels.tsv   — Sample \t Cluster (Mn/Mf/Peninsular) for the 501 ∩
#   $OUT/Exclude_for_Mf.tsv   — plink --remove format (all NON-Mf samples in the 501)
#   $OUT/Exclude_for_Mn.tsv
#   $OUT/Exclude_for_Pen.tsv
for f in keep_501.txt cluster_labels.tsv Exclude_for_Mf.tsv Exclude_for_Mn.tsv Exclude_for_Pen.tsv; do
  [[ -s "$OUT/$f" ]] || { echo "ERROR: expected $OUT/$f missing after 12b" >&2; exit 1; }
done
for c in Mf Mn Pen; do
  echo "  $c cluster kept = $(($(wc -l < "$OUT/keep_501.txt") - $(wc -l < "$OUT/Exclude_for_$c.tsv")))"
done

# ---- (2) subset regen cleaned to the 501 ∩ ----------------------------------
echo "--- plink --keep restrict to 501 ∩ ---"
plink --bfile "$PLINK/cleaned" \
      --keep "$OUT/keep_501.txt" \
      --make-bed \
      --allow-no-sex \
      --out "$OUT/cleaned_501"
wc -l "$OUT/cleaned_501.fam" "$OUT/cleaned_501.bim"

# ---- (3) per-cluster --remove Exclude → --maf 0.05 ---------------------------
for C in Mf Mn Pen; do
  echo "--- plink per-cluster $C: remove non-$C, --maf 0.05 ---"
  plink --bfile "$OUT/cleaned_501" \
        --remove "$OUT/Exclude_for_$C.tsv" \
        --recode \
        --allow-no-sex \
        --out "$OUT/$C"
  plink --file "$OUT/$C" \
        --maf 0.05 \
        --make-bed \
        --allow-no-sex \
        --out "$OUT/$C"
  wc -l "$OUT/$C.fam" "$OUT/$C.bim"
done

# ---- (4) merge per-cluster .bed → clusters_combined -------------------------
echo "--- plink --merge-list → clusters_combined ---"
ls "$OUT"/{Mf,Mn,Pen}.bed | sed 's/\.bed$//' > "$OUT/cluster_files.txt"
cat "$OUT/cluster_files.txt"
plink --merge-list "$OUT/cluster_files.txt" \
      --make-bed \
      --allow-no-sex \
      --out "$OUT/clusters_combined"
wc -l "$OUT/clusters_combined.fam" "$OUT/clusters_combined.bim"

# ---- (5) restrict to GIS-non-missing (mirrors L1266) + sed */0 → *_fixed.ped ----
echo "--- restrict to GIS-non-missing samples (via complete_gis_na_omit.fam) + sed */0 ---"
awk '{print $1"\t"$2}' data/raw/tess3r/complete_gis_na_omit.fam > "$OUT/keep_gis.txt"
plink --bfile "$OUT/clusters_combined" \
      --recode \
      --keep "$OUT/keep_gis.txt" \
      --tab \
      --allow-no-sex \
      --out "$OUT/clusters_combined"
sed 's/\*/0/g' "$OUT/clusters_combined.ped" > "$OUT/clusters_combined_fixed.ped"

echo "--- DONE ---"
echo "clusters_combined.fam           = $(wc -l < $OUT/clusters_combined.fam)  (poster 532)"
echo "clusters_combined.map           = $(wc -l < $OUT/clusters_combined.map)  (poster ~75,912)"
echo "clusters_combined_fixed.ped rows= $(wc -l < $OUT/clusters_combined_fixed.ped)"
echo "clusters_combined_fixed.ped cols= $(awk 'NR==1{print NF; exit}' $OUT/clusters_combined_fixed.ped)  (poster 151,830 = 6 + 2×~75,912)"
