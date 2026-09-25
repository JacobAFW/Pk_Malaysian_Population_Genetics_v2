# Phase 3 — deriving the support floor on this cohort, both arms

**Scope (reviewer's phase-2 spec).** Derive the per-cluster support floor by permutation with an
explicit FDR target, on this cohort, for both arms; then run the stability sweep around the
chosen N. Read the cached calls; do not re-run detection.

**Done-condition.** Four non-empty TSVs; chosen N, binding cluster, FDR at chosen N and
expected-null-windows all traceable; if Peninsular retains zero windows at the chosen N, report
it as a finding and tag NEEDS-JACOB — do not tune, do not switch to per-cluster floors. — **Met.**

---

## What was run

```
scripts/introgression_module/introgression_floor_derivation.R
  --window-size 10000 --min-snps 5 --min-samples-per-window 2
  --permutations 1000 --seed 20260925
  --max-n 60 --target-false-windows 1 --fdr-target 0.05
  --clusters <arm cluster table>   <arm's three cached call tables>

scripts/introgression_module/introgression_floor_sweep.R
  --n-values "1,3,5,6,10,16,20,25,30,35,38,39,40,41,42,44,50"
  --gff data/gadi/gff/... --fai data/gadi/fasta/... --gene-family-filters "SICA,KIR"
```

**Parameter sources, all stated.** `permutations 1000`, `target-false-windows 1`,
`fdr-target 0.05`, `min-samples-per-window 2`, window 10 kb and `min-snps 5` are the module
defaults documented in the script header and agnostic spec §9.5. **Seed 20260925** is this run's
own, declared in `config/malay_cohort.yaml` and used identically for both arms. The one value I
had to choose is **`--max-n 60`**: the module default is 25, which is *below* the answer on this
cohort and would have truncated the search — the spec records N = 44 for the 558-basis Malay
benchmark, so any bound under that returns the bound instead of the derivation (the script does
warn when no N within the bound meets the target, so it would not have been silent). 60 is a
search ceiling, not a threshold; the chosen N (41/40) sits well inside it.

**Arm-specific cluster tables.** The floor's null scales with each cluster's true n, so the
de-clonalized arm was given `inputs/admix_clusters_declonal.tsv` (484: Mf 359 · Mn 97 ·
Peninsular 28) rather than reusing the full table.

## Outputs

`results/<run>/floor/floor_derivation_{full,declonal}.tsv` and
`floor_stability_sweep_{full,declonal}.tsv` (+ `_sweep_{full,declonal}/` work dirs).

## Key numbers

**Eligibility universe** (recomputed from the genotype table, not assumed):
**1,363 windows**; per-sample eligible median 1,362 (min 1,222). Full arm **53,615** raw calls →
**46,954** distinct (sample, window) over 501 samples; de-clonalized 50,970 → 44,662 over 484. The script's
own consistency check passed on both arms: every observed call lies inside the eligible set.
Source: `floor_derivation_*.tsv` columns `universe_windows`,
`eligible_windows_per_sample_median`.

**Null support distribution and the FDR-safe floor per cluster** — from `floor_derivation_*.tsv`:

| arm | cluster | n | null mean | p95 | p99 | max | FDR-safe N | expected null windows | observed | FDR |
|---|---|---|---|---|---|---|---|---|---|---|
| full | Mf | 366 | 23.80 | 32 | 35 | 50 | **41** | 0.795 | 264 | 0.0030 |
| full | Mn | 105 | 7.68 | 12 | 14 | 24 | 18 | 0.672 | 189 | 0.0036 |
| full | Peninsular | 30 | 3.09 | 6 | 7 | 14 | 10 | 0.398 | 89 | 0.0045 |
| declonal | Mf | 359 | 22.89 | 31 | 34 | 48 | **40** | 0.674 | 258 | 0.0026 |
| declonal | Mn | 97 | 6.80 | 11 | 13 | 21 | 17 | 0.441 | 174 | 0.0025 |
| declonal | Peninsular | 28 | 3.19 | 6 | 7 | 14 | 10 | 0.547 | 109 | 0.0050 |

**Chosen global floor: N = 41 (full), N = 40 (de-clonalized); binding cluster Mf in both.**
Worst-case expected null windows 0.795 / 0.674 (target < 1); worst-case FDR 0.0030 / 0.0026
(target < 0.05).

**At the chosen floor** (`is_chosen` rows): Mf 264 / 258 windows, Mn 22 / 17, **Peninsular 0 / 0**.

**Stability sweep** — windows surviving the *full* filter chain (support floor, then gene-family
and hypervariable masks), from `floor_stability_sweep_*.tsv`:

| N | full: Mf | Mn | Pen | all | declonal: Mf | Mn | Pen | all |
|---|---|---|---|---|---|---|---|---|
| 5 | 252 | 123 | 72 | 447 | 275 | 104 | 79 | 458 |
| 10 | 354 | 140 | 40 | **534** | 360 | 125 | 45 | **530** |
| 16 | 350 | 104 | 4 | 458 | 360 | 86 | 2 | 448 |
| 20 | 336 | 80 | 1 | 417 | 333 | 61 | **0** | 394 |
| 25 | 306 | 49 | **0** | 355 | 306 | 40 | 0 | 346 |
| 35 | 271 | 25 | 0 | 296 | 267 | 17 | 0 | 284 |
| **40** | 255 | 15 | 0 | 270 | **244** | **10** | 0 | **254** |
| **41** | **249** | **14** | 0 | **263** | 242 | 7 | 0 | 249 |
| 44 | 238 | 10 | 0 | 248 | 231 | 2 | 0 | 233 |
| 50 | 217 | 2 | 0 | 219 | 203 | 0 | 0 | 203 |

Two counts differ by design and should not be confused: the derivation's `observed_windows`
(Mf 264 at N=41) counts windows passing the **support floor only**; the sweep's `n_windows`
(Mf 249) counts windows passing the **whole filter chain**, including the gene-family and
hypervariable masks that run *after* the floor. Both are in the tables named above.

**Cost.** Derivation 10.0 s / 3.3 GB (full) and 9.6 s (declonal); sweeps 30.8 s and 30.3 s at
0.32 GB. Nothing near the 2 h limit.

## Findings

**1. The Indo §9.5 pattern reproduces on this cohort, and it is a finding, not a setting.**
The FDR-safe floor scales with cluster size — Mf 41, Mn 18, Peninsular 10 (full arm) — so a
single global floor is set by the largest cluster and is four times what the smallest cluster's
own FDR argument requires. At N = 41 the Peninsular cluster retains **zero** windows. It does
not die at the floor: on the sweep grid it is already at zero from N = 25 in the full arm and
N = 20 de-clonalized, and the derivation table (support floor only, no grid) puts the onset
earlier still at **N = 23 / N = 18** — all far below 41. Per mission §5 this is reported, not tuned around, and
per-cluster floors (Mf 41 / Mn 18 / Pen 10) were **not** adopted — the module docs record that
option as explicitly rejected for reintroducing the size dependence the floor exists to remove.

**2. The mission's headline pair survives, but the Mn side is thin.** Mf↔Mn is the headline
(mission §1), and at each arm's *own* derived floor it stands: **Mf 249 windows at N = 41 (full)
and Mf 244 at N = 40 (de-clonalized)** after the full filter chain, with Mn holding **14** and
**10** respectively. So the headline is supportable, but its Mn half rests on low-double-digit
window counts, and the de-clonalized Mn column reaches zero by N = 50. This should temper how strongly a
"Mf↔Mn intermixing" claim is worded at the window level; the sample-level result from phase 2 is
the robust part.

**3. Why Peninsular cannot be rescued by the floor.** The detector calls liberally — phase 2
showed *every* eligible sample is called in *every* pair. With 366 Mf samples drawing from a
1,363-window universe, chance alone puts a mean of 23.8 Mf samples on every window (null mean,
table above). A floor strong enough to beat that (41) is far above anything a 30-sample cluster
can supply. This is the same structural argument the module's spec makes for Aceh, arrived at
independently on this cohort's own numbers.

**4. The non-monotonicity reproduces.** Total windows rise to 534 (N = 10) then fall — raising
the floor strips one cluster's weak support for a window, and a window called in two clusters
becomes called in one and *survives* the hypervariable filter that runs after. "Raise the floor
to be safer" is not a one-way lever here either.

**5. De-clonalization barely moves the floor.** 41 → 40, and the per-cluster FDR-safe values
shift by at most one (Mn 18 → 17, Peninsular 10 → 10). Removing 17 redundant genotypes changes
the null's scale slightly but not the structure of the problem.

## NEEDS-JACOB

**The FDR argument and any Peninsular window-level claim are incompatible on this cohort.**
At the FDR-safe global floor (41 full / 40 de-clonalized, binding cluster Mf), Peninsular retains
zero windows in both arms; its own cluster-scaled FDR-safe floor would be 10, where it retains
89 / 109 windows. The options, in the order the module's spec ranks them, are: (1) accept the
floor and retire Peninsular window-level claims, keeping the cross-cluster Mf/Mn result;
(2) test Peninsular at its own scale with its own null — the analogue of the focal test, which
mission §5 puts out of scope for this run; (3) reduce the detection call rate first and
re-derive; (4) per-cluster floors, which the module rejects. **Not tuned around, and no option
taken.** This is a science call.

## Proposed next (phase 4) and alternatives

Aggregate both arms at their derived floors (41 / 40) through `introgression_aggregate.R` —
the four-filter chain — writing per-pair window tables and per-sample tables, with the Mf↔Mn
headline called out and the clonality effect quantified on the *filtered* sets. Alternatives
considered: (a) aggregate at a floor that keeps Peninsular alive (e.g. 10) — rejected, that is
exactly the tuning mission §5 forbids; the sweep already reports what every floor yields, which
is the honest way to show it. (b) Aggregate only the full arm and treat de-clonalized as a
sensitivity check — rejected, the mission requires both arms reported, not one with a footnote.

---

## Revision note (phase-3 review returned REVISE)

Four corrections to this report's text. **No output changed and nothing was re-run** — the four
TSVs are unaltered and were independently verified correct; these were transcription errors on
my side.

1. **Full-arm call counts were wrong.** I wrote "50,970 raw calls → 47,079 distinct" for the full
   arm. 50,970 is the *de-clonalized* raw count, pasted into the full-arm slot, and 47,079 was not
   reproducible from any dedup key. The full arm is **53,615 raw → 46,954 distinct (sample,
   window)**; the de-clonalized figures (50,970 → 44,662) were right. Both recomputed directly
   from `calls/<arm>/*.tsv`.
2. **Finding 2 quoted the de-clonalized arm at the wrong floor.** I gave "Mf 242 / Mn 7", which is
   the sweep's N = 41 row; the de-clonalized arm's *own* chosen floor is N = 40, where the values
   are **Mf 244 / Mn 10**. This slightly strengthens the headline rather than weakening it, and
   "single-digit" was therefore wrong — corrected to low-double-digit.
3. **"Silently returns the bound" overstated the `--max-n` risk.** The script warns when no N
   within the bound meets the target, so a truncated search would have been visible. The reason
   for raising the bound stands; the characterisation was too strong.
4. **Added the earlier zero-onset for Peninsular** from the derivation table (support floor only):
   N = 23 full / N = 18 de-clonalized, earlier than the N = 25 / N = 20 the sweep grid shows, and
   more conservative for the finding.

The reviewer also verified the full four-filter chain end-to-end against the sweep's own
`aggregate.log` (53,615 → 53,584 → 26,275 → 25,520 → 24,658 calls / 263 windows / 471 samples at
N = 41), which matches the sweep's N = 41 "all" row and confirms the derivation-vs-sweep count
distinction described above.
