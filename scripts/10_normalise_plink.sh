#!/usr/bin/env bash
# Stage B-2: bcftools norm + PLINK pipeline (lifted from scripts/gadi/Analyses.Rmd).
# Inputs : data/gadi/Consensus_SNPs_subset_no_MOI.vcf.gz   (755 samples)
#          data/gadi/fasta/strain_A1_H.1.Icor.fasta(.fai)
#          data/gadi/regions_to_mask.list
# Outputs: data/processed/plink/Pk.{bed,bim,fam}            (post-dedup)
#          data/processed/plink/cleaned.{bed,bim,fam}       (target n=558)
# Run    : source vvg-box-pixi/bin/activate; bash scripts/10_normalise_plink.sh
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$(pwd)"
VCF=data/gadi/Consensus_SNPs_subset_no_MOI.vcf.gz
FASTA=data/gadi/fasta/strain_A1_H.1.Icor.fasta
MASK_SRC=data/gadi/regions_to_mask.list
OUT=data/processed/plink
mkdir -p "$OUT"

command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH; source vvg-box-pixi/bin/activate" >&2; exit 1; }
command -v plink    >/dev/null || { echo "ERROR: plink not on PATH" >&2; exit 1; }
bcftools --version | head -1
plink --version

# ---- chrom rename map (ordered_PKNH_NN_v2 → NN) ------------------------------
echo "--- chrom rename map ---"
RENAME=$OUT/chrom_rename.tsv
awk 'BEGIN{for (i=1;i<=14;i++) printf "ordered_PKNH_%02d_v2\t%d\n", i, i}' > "$RENAME"
cat "$RENAME"

# ---- regions_to_mask in plink range format (CHR BP1 BP2) ---------------------
# regions_to_mask.list has form "ordered_PKNH_NN_v2:start-end"; convert.
MASK=$OUT/regions_to_mask.txt
awk 'BEGIN{FS="[:-]"} {c=$1; sub(/ordered_PKNH_/,"",c); sub(/_v2/,"",c); printf "%d\t%d\t%d\tmask_%d\n", c+0, $2, $3, NR}' \
  "$MASK_SRC" > "$MASK"
echo "regions_to_mask.txt rows = $(wc -l < $MASK) (source $(wc -l < $MASK_SRC))"

# ---- bcftools norm (slowest step) --------------------------------------------
BCF=$OUT/Consensus_SNPs_subset_no_MOI.bcf
if [[ -s "$BCF" ]]; then
  echo "--- norm: $BCF already exists (skipping); rm to redo ---"
else
  echo "--- bcftools norm (split multiallelic + check-ref + rename chrs + annotate ID) ---"
  bcftools view "$VCF" \
    | sed 's/AD,Number=R/AD,Number=./' \
    | bcftools norm -m-any \
    | bcftools norm --check-ref w -f "$FASTA" \
    | bcftools annotate --rename-chrs "$RENAME" \
    | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT:%ID' \
    > "$BCF"
fi
ls -la "$BCF"

# ---- sample_order (preserve VCF sample order in PLINK fam) -------------------
SORDER=$OUT/sample_order.txt
bcftools query -l "$BCF" | awk '{print $1"\t"$1}' > "$SORDER"
echo "n_samples = $(wc -l < $SORDER)"   # expect 755

# ---- PLINK make-bed -> Pk ----------------------------------------------------
# Skip if Pk.{bed,bim,fam} already exist with the expected 755 (pre-metadata) rows.
if [[ -s "$OUT/Pk.bed" && -s "$OUT/Pk.bim" && -s "$OUT/Pk.fam" && \
      "$(wc -l < $OUT/Pk.fam)" == "755" ]]; then
  echo "--- plink make-bed: Pk.* already in pre-metadata state (755 fam rows); skipping ---"
else
  echo "--- plink make-bed → $OUT/Pk ---"
  rm -f "$OUT"/Pk-temporary.* "$OUT"/Pk~.*
  plink --bcf "$BCF" \
        --keep-allele-order \
        --vcf-idspace-to _ \
        --double-id \
        --allow-extra-chr 0 \
        --make-bed \
        --indiv-sort file "$SORDER" \
        --exclude range "$MASK" \
        --allow-no-sex \
        --out "$OUT/Pk"
fi
wc -l "$OUT/Pk.fam" "$OUT/Pk.bim"
echo "(checkpoint: expect 755 fam rows; bim count = retained variants after mask)"

# ---- attach metadata to Pk.fam + identify duplicates -------------------------
echo "--- attach metadata + dedup ---"
Rscript scripts/10b_attach_metadata.R

# remove duplicates → overwrites Pk.{bed,bim,fam} in $OUT (per Analyses.Rmd L586)
plink --bfile "$OUT/Pk" --recode --make-bed --remove "$OUT/Pk.dups" --allow-no-sex --out "$OUT/Pk"
echo "Pk.fam after dedup: $(wc -l < $OUT/Pk.fam)"

# ---- missingness + MAF filters → cleaned -------------------------------------
echo "--- plink filters → cleaned ---"
plink --bfile "$OUT/Pk" --recode --allow-no-sex --out "$OUT/Pk"
plink --file "$OUT/Pk" --geno 0.20 --make-bed --allow-no-sex --out "$OUT/cleaned"
plink --bfile "$OUT/cleaned" --recode --allow-no-sex --out "$OUT/cleaned"
plink --file "$OUT/cleaned" --mind 0.20 --make-bed --allow-no-sex --out "$OUT/cleaned"
plink --bfile "$OUT/cleaned" --recode --allow-no-sex --out "$OUT/cleaned"
plink --file "$OUT/cleaned" --maf 0.01 --make-bed --allow-no-sex --out "$OUT/cleaned"

echo "--- DONE ---"
echo "cleaned.fam rows = $(wc -l < $OUT/cleaned.fam)   (target = 558)"
echo "cleaned.bim rows = $(wc -l < $OUT/cleaned.bim)"
