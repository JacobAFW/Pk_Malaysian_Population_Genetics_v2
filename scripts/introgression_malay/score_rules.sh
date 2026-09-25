#!/usr/bin/env bash
# score_rules.sh — floor derivation + aggregation for every rule x scenario,
# under an explicitly stated floor policy. Reads what run_rule_comparison.sh
# produced; writes only under results/<run>/rules/<arm>/.
#
# FLOOR POLICY (stated, not implicit — a rule scored under a floor nobody chose
# is a non-result):
#   primary  "common"   — every rule is aggregated at the floor DERIVED IN
#                         PHASE 3 UNDER `absolute` (41 full / 40 declonal).
#                         The rule is then the only variable that changes, so
#                         differences are attributable to the rule.
#   cross-check "floor6" — every rule is aggregated at the PUBLISHED analysis's
#                         per-cluster floor (n > 5, i.e. 6), the configuration the
#                         Indo 558-basis scoreboard used, so our rule ranking is
#                         directly comparable to it. `full` scenario only.
#   secondary "own"     — every rule is aggregated at a floor re-derived from
#                         its OWN calls by the same permutation/FDR procedure
#                         (1,000 reps, same seed, same targets). This stops a
#                         rule with a different call rate from being judged
#                         against a floor calibrated for another rule.
# Both are reported. Neither is tuned per rule beyond these two stated rules.
#
#   bash scripts/introgression_malay/score_rules.sh [config] [arm]
set -euo pipefail

CFG="${1:-config/malay_cohort.yaml}"
ARM="${2:-full}"
cfgget() { python3 -c "
import yaml
d=yaml.safe_load(open('$CFG'))
for k in '$1'.split('.'): d=d[k]
print(d)
"; }

RES=$(cfgget outputs.results_dir)
GT=$(cfgget inputs.genotype)
WIN=$(cfgget introgression.window_size_bp)
MINSNP=$(cfgget introgression.min_snps_per_window)
MINSAMP=$(cfgget introgression.min_samples_per_window)
NPERM=$(cfgget floor_derivation.replicates)
FSEED=$(cfgget floor_derivation.seed)
FDR=$(cfgget floor_derivation.fdr_target)
TFW=$(cfgget floor_derivation.max_expected_null_windows)

case "$ARM" in
  full)     CL="$RES/inputs/admix_clusters.tsv";          COMMON_N=41;;
  declonal) CL="$RES/inputs/admix_clusters_declonal.tsv"; COMMON_N=40;;
esac

OUT="$RES/rules/$ARM"
RULES=(absolute relative distance distance_adaptive)
SCENARIOS=(full drop10 drop20)
FLOORTAB="$OUT/rule_floors.tsv"
[ -f "$FLOORTAB" ] || printf 'arm\trule\tderived_N\tbinding_cluster\n' > "$FLOORTAB"

for RULE in "${RULES[@]}"; do
  # --- per-rule derived floor, from that rule's own full-scenario calls ----
  FD="$OUT/$RULE/floor_derivation.tsv"
  if [ ! -s "$FD" ]; then
    Rscript scripts/introgression_module/introgression_floor_derivation.R \
      --genotype-table "$GT" --clusters "$CL" \
      --window-size "$WIN" --min-snps "$MINSNP" --min-samples-per-window "$MINSAMP" \
      --permutations "$NPERM" --seed "$FSEED" \
      --max-n 60 --target-false-windows "$TFW" --fdr-target "$FDR" \
      --out "$FD" "$OUT/$RULE/full"/*.tsv > "$OUT/$RULE/floor_derivation.log" 2>&1
    OWN_N=$(awk -F'\t' 'NR==2{print $11}' "$FD")
    BIND=$(awk -F'\t' 'NR==2{print $13}' "$FD")
    printf '%s\t%s\t%s\t%s\n' "$ARM" "$RULE" "$OWN_N" "$BIND" >> "$FLOORTAB"
    echo "[score_rules] $RULE: own derived floor N=$OWN_N (binding $BIND)"
  fi
  OWN_N=$(awk -F'\t' -v r="$RULE" '$2==r{print $3; exit}' "$FLOORTAB")

  # --- aggregate every scenario under both floor policies -----------------
  for SCEN in "${SCENARIOS[@]}"; do
    for POLICY in common own floor6; do
      # floor6 is the published analysis's per-cluster floor (n > 5). Only the
      # `full` scenario is aggregated there — it exists for the cross-cohort
      # comparison against the Indo 558-basis scoreboard, not for stability.
      [ "$POLICY" = "floor6" ] && [ "$SCEN" != "full" ] && continue
      case "$POLICY" in common) N=$COMMON_N;; own) N=$OWN_N;; floor6) N=6;; esac
      DEST="$OUT/$RULE/agg_${POLICY}/$SCEN"
      [ -s "$DEST/introgressed_windows_filtered.tsv" ] && continue
      mkdir -p "$DEST"
      Rscript scripts/introgression_module/introgression_aggregate.R \
        --clusters "$CL" --metadata "$RES/inputs/samples.tsv" \
        --contig-map "$RES/inputs/contig_map.tsv" \
        --fai data/gadi/fasta/strain_A1_H.1.Icor.fasta.fai \
        --gff data/gadi/gff/strain_A1_H.1.Icor.gff3 \
        --gene-family-filters "SICA,KIR" \
        --window-size "$WIN" --min-samples-per-window "$MINSAMP" \
        --per-cluster-min-pct 0 --per-cluster-min-samples "$N" \
        --out-dir "$DEST" "$OUT/$RULE/$SCEN"/*.tsv > "$DEST/aggregate.log" 2>&1
      echo "[score_rules] aggregated $RULE/$SCEN policy=$POLICY floor=$N"
    done
  done
done
echo "[score_rules] done — floors: $FLOORTAB"
