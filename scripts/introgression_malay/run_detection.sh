#!/usr/bin/env bash
# run_detection.sh — per-pair introgression detection for every pair x arm.
#
# Thin driver over the lifted scripts/introgression_module/introgression_pair.R.
# All cohort facts come from config/malay_cohort.yaml; nothing is hardcoded here
# beyond the config path default. Existing call tables are NOT recomputed.
#
#   bash scripts/introgression_malay/run_detection.sh [config] [rule]
#
# Writes only under the config's results_dir. Logs per-run wall time and peak
# RSS to results/<run>/calls/detection_log.tsv.
set -euo pipefail

CFG="${1:-config/malay_cohort.yaml}"
RULE="${2:-absolute}"

yq() { # minimal scalar reader: yq <dotted.key> — avoids a new dependency
  python3 -c "
import sys, yaml
d = yaml.safe_load(open('$CFG'))
for k in '$1'.split('.'):
    d = d[k]
print(d)
"
}

RES=$(yq outputs.results_dir)
GT=$(yq inputs.genotype)
WIN=$(yq introgression.window_size_bp)
MINSNP=$(yq introgression.min_snps_per_window)
LVL_O=$(yq introgression.contour_level_other)
LVL_W=$(yq introgression.contour_level_own)

PAIRS=(Mf__Mn Mf__Peninsular Mn__Peninsular)
LOG="$RES/calls/detection_log.tsv"
mkdir -p "$RES/calls/full" "$RES/calls/declonal"
[ -f "$LOG" ] || printf 'arm\tpair\trule\tstatus\twall_s\tpeak_rss_gb\tout\n' > "$LOG"

for ARM in full declonal; do
  case "$ARM" in
    full)     EXCL="$RES/inputs/exclude_full.txt" ;;
    declonal) EXCL="$RES/inputs/exclude_unique.txt" ;;
  esac
  for PAIR in "${PAIRS[@]}"; do
    OUT="$RES/calls/$ARM/$PAIR.tsv"
    if [ -s "$OUT" ]; then
      echo "[run_detection] $ARM/$PAIR: already on disk, reusing ($(wc -l < "$OUT") lines)"
      printf '%s\t%s\t%s\treused\tNA\tNA\t%s\n' "$ARM" "$PAIR" "$RULE" "$OUT" >> "$LOG"
      continue
    fi
    echo "[run_detection] $ARM/$PAIR: running"
    TIMEF=$(mktemp)
    /usr/bin/time -l Rscript scripts/introgression_module/introgression_pair.R \
      --genotype-table "$GT" \
      --clusters "$RES/inputs/admix_clusters.tsv" \
      --exclude "$EXCL" \
      --pair "$PAIR" \
      --window-size "$WIN" \
      --min-snps "$MINSNP" \
      --detection-rule "$RULE" \
      --contour-level-other "$LVL_O" \
      --contour-level-own "$LVL_W" \
      --out "$OUT" 2> "$TIMEF" || { cat "$TIMEF"; exit 1; }
    grep -vE '^ *[0-9]+ +[a-z]' "$TIMEF" || true
    WALL=$(awk '/real/{print $1}' "$TIMEF" | tail -1)
    RSS=$(awk '/maximum resident set size/{printf "%.1f", $1/1073741824}' "$TIMEF" | tail -1)
    printf '%s\t%s\t%s\tok\t%s\t%s\t%s\n' "$ARM" "$PAIR" "$RULE" "$WALL" "$RSS" "$OUT" >> "$LOG"
    rm -f "$TIMEF"
  done
done

echo
echo "[run_detection] done — log: $LOG"
