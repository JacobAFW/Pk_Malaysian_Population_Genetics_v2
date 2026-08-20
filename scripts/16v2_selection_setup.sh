#!/usr/bin/env bash
# Stage B-8 v2 setup: per-cluster VCFs on the FULL regen cohort (559) with
# regen ADMIXTURE labels (drops the 501 ∩ F-anchor for selection — cohort
# size matters for iHS power).
set -euo pipefail

cd "$(dirname "$0")/.."
PLINK=data/processed/plink
NO_MOI=data/gadi/Consensus_SNPs_subset_no_MOI.vcf.gz
OUT=data/processed/selection_v2
mkdir -p "$OUT"

command -v plink >/dev/null    || { echo "ERROR: plink not on PATH" >&2; exit 1; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH" >&2; exit 1; }

# ---- (1) regen cluster labels + per-cluster exclude files -------------------
echo "--- derive regen cluster labels from regen cleaned.3.Q (559) ---"
Rscript scripts/16v2_regen_labels.R

for c in Mf Mn Pen; do
  [[ -s "$OUT/Exclude_for_$c.tsv" ]] || { echo "ERROR: expected $OUT/Exclude_for_$c.tsv missing" >&2; exit 1; }
done

# ---- (2) per-cluster PLINK subset + --maf 0.05 (regen-cohort BIM per cluster) ----
for C in Mf Mn Pen; do
  echo "--- plink per-cluster $C: --remove non-$C + --maf 0.05 (from regen cleaned.559) ---"
  plink --bfile "$PLINK/cleaned" \
        --remove "$OUT/Exclude_for_$C.tsv" \
        --recode \
        --allow-no-sex \
        --out "$OUT/${C}_maf_step1"
  plink --file "$OUT/${C}_maf_step1" \
        --maf 0.05 \
        --make-bed \
        --allow-no-sex \
        --out "$OUT/$C"
  rm -f "$OUT/${C}_maf_step1".{ped,map,log,nosex,irem}
  echo "  $C.fam n = $(wc -l < $OUT/$C.fam) | $C.bim n = $(wc -l < $OUT/$C.bim)"
done

# ---- (3) per-cluster VCFs + per-chrom split (same as v1, using v2 labels) --
CONTIGS=(ordered_PKNH_01_v2 ordered_PKNH_02_v2 ordered_PKNH_03_v2 ordered_PKNH_04_v2 \
         ordered_PKNH_05_v2 ordered_PKNH_06_v2 ordered_PKNH_07_v2 ordered_PKNH_08_v2 \
         ordered_PKNH_09_v2 ordered_PKNH_10_v2 ordered_PKNH_11_v2 ordered_PKNH_12_v2 \
         ordered_PKNH_13_v2 ordered_PKNH_14_v2)

for CLUSTER in Mf Mn Pen; do
  CDIR=$OUT/${CLUSTER}_vcf
  mkdir -p "$CDIR"

  # snp_filter.txt from per-cluster .bim (chrom \t pos). BIM uses numeric chroms
  # (1..14) — need to prefix "ordered_PKNH_NN_v2" for VCF filtering.
  awk 'BEGIN{OFS="\t"} {printf "ordered_PKNH_%02d_v2\t%s\n", $1, $4}' \
      "$OUT/$CLUSTER.bim" > "$CDIR/snp_filter.txt"
  echo "--- $CLUSTER: snp_filter lines = $(wc -l < $CDIR/snp_filter.txt) ---"

  # Filter VCF: only samples in this cluster (using v2 Exclude), only snps in filter
  if [[ ! -s "$CDIR/selection.vcf.gz" ]]; then
    echo "  bcftools view + missing2ref → $CDIR/selection.vcf.gz"
    bcftools view --force-samples \
             -R "$CDIR/snp_filter.txt" \
             -S "^$OUT/Exclude_for_${CLUSTER}_bcftools.txt" \
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
  echo "  per-chrom VCFs = $(ls $CDIR/ordered_PKNH_*.vcf.gz 2>/dev/null | wc -l)"
done

echo "--- DONE ---"
