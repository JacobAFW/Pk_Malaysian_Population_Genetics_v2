#!/usr/bin/env bash
# run_rule_comparison.sh — detection under every rule x perturbation scenario.
#
# Phase 6 of the run. Produces the raw material the scorer reads; it does not
# score anything and it does not touch the headline outputs under
# results/<run>/{calls,aggregate,floor}.
#
# Scenarios are sample-DROP perturbations expressed as exclusion lists, so the
# same `--exclude` machinery the de-clonalized arm uses drives them, and
# cluster DEFINITIONS always come from the full table (module convention).
#
#   bash scripts/introgression_malay/run_rule_comparison.sh [config] [arm]
#
# Existing outputs are reused, never recomputed.
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
LVL_O=$(cfgget introgression.contour_level_other)
LVL_W=$(cfgget introgression.contour_level_own)
MARGIN=$(cfgget introgression.distance_margin)
ADQ=$(cfgget introgression.distance_adaptive_quantile)
SEED=$(cfgget benchmark.seed)

PAIRS=(Mf__Mn Mf__Peninsular Mn__Peninsular)
RULES=(absolute relative distance distance_adaptive)
SCENARIOS=(full drop10 drop20)

OUT="$RES/rules/$ARM"
mkdir -p "$OUT"
LOG="$OUT/rule_detection_log.tsv"
[ -f "$LOG" ] || printf 'arm\trule\tscenario\tpair\tstatus\twall_s\n' > "$LOG"

# Base cluster table for the arm, and the arm's own exclusion list.
case "$ARM" in
  full)     BASE_CL="$RES/inputs/admix_clusters.tsv";          ARM_EXCL="";;
  declonal) BASE_CL="$RES/inputs/admix_clusters_declonal.tsv"; ARM_EXCL="$RES/inputs/exclude_unique.txt";;
  *) echo "unknown arm: $ARM" >&2; exit 1;;
esac

# --- perturbation exclusion lists (deterministic, one stated seed) ---------
# Drawn once per arm and reused by every rule, so rules see identical cohorts.
python3 - "$BASE_CL" "$OUT" "$SEED" "${ARM_EXCL:-}" <<'PY'
import csv, random, sys, os
base, out, seed, arm_excl = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
samples = [r["Sample"] for r in csv.DictReader(open(base), delimiter="\t")]
samples.sort()
carried = []
if arm_excl and os.path.exists(arm_excl):
    carried = [l.strip() for l in open(arm_excl) if l.strip()]
for frac, name in ((0.10, "drop10"), (0.20, "drop20")):
    p = os.path.join(out, f"exclude_{name}.txt")
    if os.path.exists(p):
        continue
    rng = random.Random(seed)              # same seed -> same draw for every rule
    k = int(round(frac * len(samples)))
    dropped = sorted(rng.sample(samples, k))
    # the arm's own exclusions must stay excluded on top of the perturbation
    with open(p, "w") as fh:
        fh.write("\n".join(sorted(set(dropped) | set(carried))) + "\n")
    print(f"[rule_comparison] {name}: dropped {k} of {len(samples)} (seed {seed})")
PY

for RULE in "${RULES[@]}"; do
  case "$RULE" in
    distance_adaptive) DET=distance; ADAPT=true  ;;
    *)                 DET="$RULE";  ADAPT=false ;;
  esac
  for SCEN in "${SCENARIOS[@]}"; do
    case "$SCEN" in
      full)   EXCL="${ARM_EXCL:-$RES/inputs/exclude_full.txt}";;
      drop10) EXCL="$OUT/exclude_drop10.txt";;
      drop20) EXCL="$OUT/exclude_drop20.txt";;
    esac
    mkdir -p "$OUT/$RULE/$SCEN"
    for PAIR in "${PAIRS[@]}"; do
      DEST="$OUT/$RULE/$SCEN/$PAIR.tsv"
      # The headline absolute/full calls already exist — link, never recompute.
      SRC="$RES/calls/$ARM/$PAIR.tsv"
      if [ "$RULE" = "absolute" ] && [ "$SCEN" = "full" ] && [ -s "$SRC" ] && [ ! -s "$DEST" ]; then
        cp "$SRC" "$DEST"
        printf '%s\t%s\t%s\t%s\treused_headline\tNA\n' "$ARM" "$RULE" "$SCEN" "$PAIR" >> "$LOG"
        continue
      fi
      if [ -s "$DEST" ]; then
        printf '%s\t%s\t%s\t%s\treused\tNA\n' "$ARM" "$RULE" "$SCEN" "$PAIR" >> "$LOG"
        continue
      fi
      T0=$(python3 -c 'import time; print(time.time())')
      Rscript scripts/introgression_module/introgression_pair.R \
        --genotype-table "$GT" \
        --clusters "$BASE_CL" \
        --exclude "$EXCL" \
        --pair "$PAIR" \
        --window-size "$WIN" \
        --min-snps "$MINSNP" \
        --detection-rule "$DET" \
        --contour-level-other "$LVL_O" \
        --contour-level-own "$LVL_W" \
        --distance-margin "$MARGIN" \
        --distance-adaptive "$ADAPT" \
        --distance-adaptive-quantile "$ADQ" \
        --out "$DEST" > "$OUT/$RULE/$SCEN/$PAIR.log" 2>&1
      T1=$(python3 -c 'import time; print(time.time())')
      printf '%s\t%s\t%s\t%s\tok\t%.1f\n' "$ARM" "$RULE" "$SCEN" "$PAIR" \
        "$(python3 -c "print($T1-$T0)")" >> "$LOG"
      echo "[rule_comparison] $ARM/$RULE/$SCEN/$PAIR done"
    done
  done
done
echo "[rule_comparison] complete — $LOG"
