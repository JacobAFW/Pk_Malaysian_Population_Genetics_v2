# RUN SUMMARY — 2026-09-25_introgression-module

Agnostic introgression module applied to the Large-Malay cohort.
Branch `auto/2026-09-25_introgression-module` (never merged, never pushed; `main` untouched).
Phases 0–7 of a maximum 8, all reviewed. Run wall-clock ≈ 1 h 30 m.

**Start here:** `results/2026-09-25_introgression-module/notes.md` — the drafting document.
Its **§7 corrections register** lists every superseded number found during the run; read it before
quoting any figure from an individual phase report, because four reports contain numbers that were
later corrected.

---

## 1. NEEDS-JACOB — decisions waiting on you

1. **Peninsular retains zero windows at the FDR-safe floor** (41 full / 40 de-clonalized, binding
   cluster Mf) in both arms; its own cluster-scaled floor would be 10, where it keeps 89 / 109
   windows. *Recommend:* accept the derived floor, retire Peninsular window-level claims, keep the
   cross-cluster Mf/Mn result; an own-scale Peninsular null is a separate run. Not tuned around —
   this is the mission §5 pre-registered incompatibility.
   → `results/<run>/floor/floor_derivation_{full,declonal}.tsv`, notes §2.
2. **Confirm the de-clonalized arm's basis.** The two clonal tables named in `MISSION.md` §3 carry
   field-collection ids with **zero** overlap against the genotype namespace — used as given they
   produce a silent no-op. The arm was instead derived from this project's own pairwise hmmIBD
   output at `fract_sites_IBD >= 0.98` with a lowest-missingness representative → 501 → 484.
   *Recommend:* approve as derived (19 of the 21 originally named clonal samples are recovered;
   fully reversible via the audit table). → `results/<run>/inputs/declonalization_audit.tsv`, notes §3.
3. **Confirm the shipped rule for the paper: `absolute` @ 5e-4 at floor 41.** *Recommend:* accept,
   with two caveats — the margin over `relative` is small and floor-contingent, and
   `distance_adaptive` wins decisively at permissive floors. → `results/<run>/benchmark/rule_scoreboard.tsv`, notes §6.
4. **A pre-existing ENA run accession appears inline in 5 `scripts/gadi/` files** (12 occurrences),
   committed before this run and untouched by it. Choose forward-redaction on `main` vs a history
   purge. Not actioned here: out of scope, and a history rewrite is a destructive git operation the
   protocol forbids.

## 2. What changed

### Phases and verdicts

| phase | what | verdict |
|---|---|---|
| 0 | plan | PROCEED — reviewer caught that the poster's filtered truth has no coordinate columns; `mn_windows.tsv` had to be sourced separately |
| 1 | lift module, assemble inputs, smoke-test | PROCEED |
| 2 | detection, 3 pairs × both arms | PROCEED |
| 3 | support-floor derivation | **REVISE → PROCEED** — two transcription errors: full-arm call counts (53,615 → 46,954, not 50,970 → 47,079) and the de-clonalized arm quoted at the wrong floor (N=40 gives Mf 244 / Mn 10, not 242 / 7) |
| 4 | aggregation at the derived floors | PROCEED |
| 5 | coordinate benchmark | **REVISE → PROCEED** — my "the floor removes noise rather than signal" mechanism was **wrong**; it is a floor × filter-4 interaction |
| 6 | detection-rule comparison | PROCEED |
| 7 | figures + notes.md | PROCEED |

Every phase report was independently recomputed by the reviewer — in later phases cell-by-cell with
an independent implementation rather than spot-checked.

### Headline results

| result | value | file |
|---|---|---|
| derived support floor | **41** (full) / **40** (de-clonalized), binding cluster Mf, FDR 0.0030 / 0.0026 | `results/<run>/floor/floor_derivation_*.tsv` |
| filtered result | 24,658 calls / **263 windows** / 471 samples (full); 22,946 / 254 / 456 (de-clonalized) | `results/<run>/aggregate/*/filter_audit.tsv` |
| per cluster | Mf 249 windows · Mn 14 · **Peninsular 0** (full) | `results/<run>/aggregate/full/windows_by_cluster.tsv` |
| headline Mf↔Mn | 222 vs 11 windows — **~20:1 one-directional**, Mf carrying Mn-like windows | `results/<run>/aggregate/full/introgressed_windows_filtered.tsv` |
| coordinate benchmark | **45.6%** (full) / 42.9% (de-clonalized) of the poster's 217 windows; Mf side **58.0%** vs a 64% published ceiling | `results/<run>/benchmark/coordinate_benchmark.tsv` |
| rule comparison | ranking **inverts with the floor**; `absolute` mid-table everywhere; recommend keeping it | `results/<run>/benchmark/rule_scoreboard.tsv` |
| clonality | sample *status* is a clean subtraction; the *call set* is not (Jaccard 0.883, 454/456 samples changed) | notes §3 |

### Outputs

- `results/2026-09-25_introgression-module/` — `inputs/` (derived tables + provenance + coordinate
  truth), `calls/` (6 per-pair tables), `floor/` (derivation + stability sweeps, both arms),
  `aggregate/` (filtered windows, per-sample summaries, filter audits, both arms), `benchmark/`
  (coordinate benchmark, context, rule scoreboard), `rules/` (4 rules × 3 scenarios × 3 pairs +
  per-rule floors and aggregations), **`notes.md`**.
- `figures/2026-09-25_introgression-module/` — 5 figures, vector PDF + 300-dpi PNG each. Existing
  `figures/` untouched.
- `scripts/introgression_module/` — the lifted module (9 files, md5-verified against source commit
  `d6a788e`; the read-only source tree was never written to).
- `scripts/introgression_malay/` — cohort drivers: `prepare_inputs.R`, `run_detection.sh`,
  `benchmark_coordinates.R`, `run_rule_comparison.sh`, `score_rules.sh`, `score_rules.py`,
  `make_figures.R`. Every result regenerates from these plus `config/malay_cohort.yaml`.

`results/` and `figures/` are gitignored; only scripts, config and `runs/` are tracked, and no
tracked file contains a sample identifier.

### Two findings worth carrying into the paper

- **The floor and filter 4 interact.** Filter 4 (multi-cluster mask) runs *after* the support
  floor, so raising the floor makes windows single-cluster and *rescues* them from the mask. At
  floor 6 that mask swallows 128 of the 217 truth windows; at floor 41, one. Any comparison across
  floors is partly mask bookkeeping — notes §5.
- **Independent replication of the published rule scoreboard.** At matched floor 6 our 501-basis
  chain reproduces the 558-basis numbers (37.3 vs 37%, 13.8 vs 14%, 5.5 vs 7%, 47.0 vs 49%). Scope
  it correctly: same code both times, so it validates the lift and shows the basis change barely
  matters — not validation against an independent implementation.

## 3. Low-confidence decisions — look at these first

The reviewer recorded these at **medium** confidence; they are the ones most likely to need revisiting.

1. **Accepting the substituted de-clonalization basis** (phase 1). The 0.98 threshold's provenance
   is *descriptive* — it appears in `scripts/14b_clonal_networks.R` as a plot title describing the
   truth table, not as an executed filter. Corroborated but not operational. → NEEDS-JACOB 2.
2. **Full-arm-only rule comparison** (phase 6). Detection, floor and benchmark are both-arm; the
   rule sweep is not. A de-clonalized rule ranking flip is unlikely but untested.
3. **The `absolute` recommendation** (phase 6). `relative` reproduces 7 more truth windows at the
   shipped floor; ~29% of that margin is mask movement, and `absolute` is more stable at drop10.
   Close call. → NEEDS-JACOB 3.

Also note: I computed F_MISS from the genotype table because no PLINK `.smiss` exists anywhere in
the project. This changes the chosen clonal representative in **8 of 15** groups versus the
documented alphabetical-only fallback. The reviewer judged it the rule faithfully applied rather
than a deviation, but it is consequential and both choices are recorded per group.

## 4. Suggested next run

1. **De-clonalized rule sweep** — closes the phase-6 coverage gap; ~15 min via the committed
   scripts (`run_rule_comparison.sh` / `score_rules.sh` / `score_rules.py` with `declonal`).
2. **Peninsular at its own scale** — an own-null test rather than filtering it through a
   cluster-wide floor. This is option 2 of the module's own list and the natural answer to
   NEEDS-JACOB 1.
3. **Introgression × environment association** (ALC vector / tree cover / elevation) and its PNTD
   replication — explicitly out of scope for this run (mission §5).
4. Optional pre-submission figure polish, from the phase-7 review: give directions and floor
   policies their own palette rather than reusing the cluster hues (purple means Peninsular in
   figs 1/3b/4 but something else in 3a and 5), and relabel fig 4a's "raw calls" as "raw window
   sets" to avoid conflating it with the call-level 0.883 in notes §3. Neither is a wrong claim.

## 5. Corrections and blocked actions

- **Correction to the phase-2 report and the phase-1/2 log entries:** the pre-existing ENA
  accession appears in **5** `scripts/gadi/` files, not 6 (verified: 12 occurrences across
  `Analyses.Rmd`, `high_quality_with_cluster_maf.Rmd`, `hmmIBD_cluster_spec.pbs`,
  `hmmIBD_no_lab_strains.pbs`, `hmmIBD_strict.pbs`). Log entries are append-only and were not
  edited.
- **Redaction (phase 2).** The reviewer's phase-01 log entry quoted two real sample ids. Since
  `REVIEW_LOG.md` is tracked in a public repo, they were replaced with markers pointing at the
  gitignored audit table *before* the entry was ever committed — the ids never entered git history.
  Recorded in the log under `## Blocked actions / redactions`.
- **Blocked action (phase 7).** The floor-6 regeneration check left a scratch copy at
  `results/<run>/rules/full/absolute/agg_floor6_prev/`. The guard denies deletes outside scratch
  dirs, so it remains; it is gitignored, harmless, and unlocks nothing. Not worked around.
- **A process failure worth knowing about (phase 2).** My pre-commit sensitive-data scan returned a
  false CLEAN twice: this shell is zsh, where unquoted `$FILES` does not word-split, and `2>/dev/null`
  hid the resulting error. Scans now iterate file-by-file. The reviewer independently confirmed the
  tracked tree carries zero sample identifiers.
- **A silent bug I introduced and caught (phase 5).** The benchmark's first version keyed the
  coordinate join with `paste0()`; R renders round doubles in scientific notation
  (`as.character(100000)` → `"1e+05"`) while the integer `START` rendered `"100000"`, so 17 truth
  coordinates were unmatchable. It under-reported reproduction as 89/217 instead of 99/217 and
  corrupted the eligibility count. Found only because the same quantity computed in Python
  disagreed. Every load-bearing number is now computed twice, in two languages.
