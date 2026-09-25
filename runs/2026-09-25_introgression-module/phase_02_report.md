# Phase 2 — cross-cluster detection, 3 pairs × both arms

**Scope (reviewer's phase-1 spec).** Run the lifted `introgression_pair.R` under `absolute`
@ 5e-4 for `Mf__Mn`, `Mf__Peninsular`, `Mn__Peninsular` on both the full (501) and de-clonalized
(484) sample sets; reuse the existing `calls/full/Mf__Mn.tsv`; confirm from the code path that
exclusion is applied *before* cluster-consensus computation.

**Done-condition.** 6 non-empty call tables, per-table counts and direction splits, a
full-vs-declonal delta per pair, wall time and peak RSS per run. — **Met.**

---

## Exclusion happens before the consensus — confirmed from the code, not assumed

In `scripts/introgression_module/introgression_pair.R`:

| line | what happens |
|---|---|
| 134–135 | the exclusion list is read and `clusters_in` is filtered: `filter(!(SAMPLE %in% drop))` |
| 141 | `pair_members` is derived **from the already-filtered** membership |
| 165 | the genotype table is subset to `want <- intersect(pair_members$SAMPLE, names(gt_wide))` |
| 180 | `dominant_alleles(gt_long, c(kx, ky))` — the consensus is computed on that subset only |

So a dropped clone contributes to neither the consensus alleles nor the density surface. This
matters because the cluster consensus is cohort-derived (agnostic spec §9.4.4): a post-hoc
sample filter would *not* have produced a de-clonalized arm, only a de-clonalized view of a
clonal analysis. The run logs confirm it at execution time — each declonal run printed
`exclusion list exclude_unique.txt: dropped 17, 484 samples remain`, then `Mf n=359, Mn n=97`
(and `Peninsular n=28`), matching `declonalization_counts.tsv` exactly.

## What was run

`bash scripts/introgression_malay/run_detection.sh config/malay_cohort.yaml absolute` — a thin
driver that reads every parameter from `config/malay_cohort.yaml`, skips any call table already
on disk, and records wall time + peak RSS per run. Parameters identical across all six runs:
window 10,000 bp, `min-snps 5`, `detection-rule absolute`, both contour levels 5e-4, 3 threads.
Arms differ only in `--exclude` (`exclude_full.txt` is empty → a guaranteed no-op).

## Outputs

`results/2026-09-25_introgression-module/calls/{full,declonal}/{Mf__Mn,Mf__Peninsular,Mn__Peninsular}.tsv`
plus `calls/detection_log.tsv` (per-run status, wall time, peak RSS).

## Key numbers

All traceable to `calls/<arm>/<pair>.tsv`; counts are independent `awk` passes over those files.

| arm | pair | calls | windows | samples | direction split |
|---|---|---|---|---|---|
| full | **Mf__Mn** | 22,371 | 1,040 | 471 | Mf_like_Mn 16,450 · Mn_like_Mf 5,921 |
| full | Mf__Peninsular | 21,031 | 917 | 396 | Mf_like_Pen 18,137 · Pen_like_Mf 2,894 |
| full | Mn__Peninsular | 10,213 | 1,028 | 135 | Mn_like_Pen 8,196 · Pen_like_Mn 2,017 |
| declonal | **Mf__Mn** | 21,672 | 1,041 | 456 | Mf_like_Mn 16,131 · Mn_like_Mf 5,541 |
| declonal | Mf__Peninsular | 20,385 | 983 | 387 | Mf_like_Pen 17,144 · Pen_like_Mf 3,241 |
| declonal | Mn__Peninsular | 8,913 | 1,041 | 125 | Mn_like_Pen 6,849 · Pen_like_Mn 2,064 |

**Clonality effect (full → de-clonalized).** Source: `calls/detection_log.tsv` + set comparisons
over the six tables.

| pair | Δ calls | Δ windows | Δ samples | window Jaccard | samples lost that are *not* dropped clones |
|---|---|---|---|---|---|
| Mf__Mn | −699 (−3.1%) | +1 (1,040→1,041) | −15 | **0.993** | 0 |
| Mf__Peninsular | −646 (−3.1%) | **+66** (917→983) | −9 | 0.917 | 0 |
| Mn__Peninsular | −1,300 (−12.7%) | +13 (1,028→1,041) | −10 | 0.930 | 0 |

**Cost.** `Mf__Mn` 38.9 s / 4.8 GB · `Mf__Peninsular` 34.3 s / 4.6 GB · `Mn__Peninsular`
13.6 s / 3.1 GB (declonal arm; full arm within 1 s of each). Five runs ≈ 2 min total, one
reused. Far inside the 2 h/step budget.

## What the numbers say so far

**1. Sample-level calls are exactly stable under de-clonalization.** In every pair, the samples
that stop being called are *precisely* the dropped clonal duplicates — 15, 9 and 10 of them, the
exact per-pair drop counts — and **zero** non-clonal samples changed status in either direction
(no losses, no gains). This independently reproduces the Stage-B diagnosis's finding that "which
samples show introgression" is the robust part of the method, and it means the de-clonalized arm
is a clean subtraction at the sample level rather than a re-analysis with a different answer.

**2. Window sets are far more stable than a random drop of similar size would suggest.**
Jaccard 0.92–0.99 against the full arm, versus 0.66 / 0.61 for the Indo benchmark's random 10% /
20% drops (agnostic spec §9.4.3). Removing 17 of 501 (3.4%) is both a smaller and a *targeted*
perturbation — it deletes redundant genotypes rather than information — so this is the expected
direction, but it is worth stating because it bounds how much of any later full-vs-declonal
difference can be attributed to cohort-shift noise.

**3. Window counts rise slightly while call counts fall.** Most visibly `Mf__Peninsular`:
−646 calls but **+66 windows**. Fewer samples cannot add calls, so the extra windows come from
the consensus and the density surface being recomputed on the de-clonalized set — removing
near-identical genotypes changes which alleles are dominant and re-shapes the contours. This is
the cohort dependence the spec documents (§9.4.4) showing up directly in our data, and it is the
reason the arm had to be built before detection rather than after. It is a method property to
report, not an error; the filtered window sets in phase 4 are where it matters.

**4. The headline pair behaves as expected.** `Mf__Mn` carries the most calls (22,371) across
1,040 windows and 471 of 501 samples — 94% of the cohort shows some Mf↔Mn signal before any
filtering. That is a very liberal raw detection rate, consistent with the spec's standing caveat
that the detector calls liberally (§9.6) and exactly why the support floor in phase 3 matters.

## Corrections and carried items

- **Correcting my phase-1 report:** I described `fract_sites_IBD >= 0.98` as this project's rule
  from `scripts/14b_clonal_networks.R`. The reviewer checked and the threshold appears there in a
  plot *title describing* the truth table, not as an executed filter — the provenance is
  descriptive, not operational. The cut still stands on the corroborating evidence (the named
  tables' own minimum 0.98672; 20 pairs above vs 3 in the preceding 0.08; 19/21 of the named
  clonal samples recovered), but it is weaker provenance than I stated, and it remains the run's
  one NEEDS-JACOB item.
- The edge case the reviewer found (a named clonal pair with only one member inside the 501, so
  no droppable group forms) is correct by construction and will be stated in the clonality notes
  in phase 7.
- Nothing was written outside `results/<run>/`; no rule other than `absolute` was run; no
  filtering was applied — phase 3 reads these raw calls.

## Proposed next (phase 3) and alternatives

Derive the per-cluster support floor on this cohort with `introgression_floor_derivation.R`
(1,000 permutations, seed 20260925, FDR target 0.05, expected-null-windows < 1), for **both
arms**, then `introgression_floor_sweep.R` for the stability sweep around the chosen N. Both read
the cached calls above; neither re-runs detection. Expect the Indo §9.5 pattern — a floor set by
the largest cluster (Mf) that Peninsular cannot reach — and report it as a finding if it appears.

Alternatives considered: (a) derive the floor on the full arm only and reuse it for the
de-clonalized arm — rejected, the null's support distribution scales with cluster size and call
rate, both of which changed, so the de-clonalized arm needs its own derivation to be honest;
(b) skip straight to per-cluster floors — rejected, the spec records that as explicitly rejected
for reintroducing size dependence, and it would pre-empt the finding rather than report it.

## Pre-commit scan — a false clean I had to correct

My first commit-guard pass on this phase reported CLEAN. It was wrong, and the reason is worth
recording because it would silently recur.

The scan built a file list into a shell variable and passed it unquoted (`grep -rnEo ... $FILES`).
**This shell is zsh, where unquoted parameter expansion does not word-split**, so grep received a
single non-existent filename; `2>/dev/null` swallowed the error and the empty output read as
"no hits". The same construction was used in the phase-1 pre-commit scan, so that clean was also
unverified.

Re-scanned every tracked file on the branch individually. Two findings:

1. **Two real sample ids (a clonal pair) were present in `REVIEW_LOG.md`**, quoted by the reviewer
   in the phase-01 entry. That file is tracked and the repo is public. Redacted before this
   commit; the redaction and its justification are recorded under `## Blocked actions / redactions`
   in the log itself.
2. **an ENA run accession appears inline in 6 files already committed before this run** (`scripts/gadi/`).
   Pre-existing, out of this run's scope, and any history rewrite is a destructive git operation
   the protocol forbids — logged NEEDS-JACOB, not acted on.

Scans from here on iterate file-by-file rather than relying on variable word-splitting.
