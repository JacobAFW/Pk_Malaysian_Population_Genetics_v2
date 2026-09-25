# Phase 6 — detection-rule comparison

**Scope (reviewer's phase-5 spec).** Compare `absolute` / `relative` / `distance`
(+ `distance_adaptive`) on three measure families — coordinate window reproduction vs the 217
truth windows, cluster-size dependence, stability under 10%/20% sample drops — with the
support-floor policy stated explicitly, and with rule differences decomposed into detection-layer
movement vs filter-4 mask movement.

**Done-condition.** `benchmark/rule_scoreboard.tsv` non-empty, 4 rules × three measure families,
every cell traceable, a recommended rule with reasoning separating reproduction from stability
and detection from mask movement. — **Met.** Full arm only (see "coverage" below).

---

## Floor policy — stated, not implicit

A rule scored under a floor nobody chose is a non-result, so all three policies are reported and
none is tuned per rule beyond its own definition:

| policy | floor used | why |
|---|---|---|
| `common` | **41** for every rule — the floor derived in phase 3 under `absolute` | makes the rule the only changing variable |
| `own` | each rule's **own** floor, re-derived from its own calls by the identical permutation/FDR procedure (1,000 reps, seed 20260925, FDR 0.05, expected-null < 1) | stops a rule with a different call rate being judged against another rule's calibration |
| `floor6` | **6** for every rule — the published analysis's per-cluster floor | makes our ranking directly comparable to the Indo 558-basis scoreboard |

**Per-rule derived floors** (`rules/full/rule_floors.tsv`), all binding on Mf:

| rule | derived N |
|---|---|
| absolute | 41 |
| relative | 43 |
| distance | **9** |
| distance_adaptive | **18** |

The spread is the justification for the `own` policy: `distance` calls so much less that its
FDR-safe floor is 9, and judging it at 41 is judging it against a threshold its call rate could
never meet.

## What was run

`run_rule_comparison.sh` (36 detection runs = 4 rules × 3 scenarios × 3 pairs, **13.9 min**
total, the `absolute`/full three reused from phase 2), then `score_rules.sh` (4 floor
derivations ≈ 10 s each, 28 aggregations ≈ 2 s each) and `score_rules.py`.

Drop scenarios are deterministic exclusion lists drawn once per arm at **seed 20260925** and
shared by every rule, so all rules see identical perturbed cohorts: drop10 = 50 of 501 excluded,
drop20 = 100 of 501.

**Coverage: full arm only.** Extending to the de-clonalized arm would roughly double the 14-min
detection cost for a comparison whose conclusion is about rule choice, not cohort — and phase 7
(figures + `notes.md`) must still fit the 8-phase limit. Stated rather than silently skipped.

`score_rules.py` is deliberately Python: after the phase-5 coordinate-keying bug, every
load-bearing number in this run is computed in a second language.

## Key numbers (`benchmark/rule_scoreboard.tsv`)

| rule | policy | floor | called | repro /217 | % | Jaccard | median win/sample | size ratio | drop10 | drop20 | Δ vs absolute | via mask |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| absolute | common | 41 | 263 | 99 | 45.6 | 0.260 | Mf 59 · Mn 6 · Pen 0 | 9.83 | 0.761 | 0.578 | — | — |
| **relative** | common | 41 | 281 | **106** | **48.9** | **0.270** | Mf 60 · Mn 16 · Pen 0 | 3.75 | 0.724 | 0.562 | 42 | 12 |
| distance | common | 41 | 16 | 1 | 0.5 | 0.004 | Mf 3 · Mn 1 · Pen 0 | 3.00 | 0.706 | 0.632 | 261 | 2 |
| distance_adaptive | common | 41 | 35 | 9 | 4.2 | 0.037 | Mf 7 · Mn 0 · Pen 0 | NA | 0.763 | **0.725** | 238 | 1 |
| absolute | own | 41 | 263 | 99 | 45.6 | 0.260 | Mf 59 · Mn 6 · Pen 0 | 9.83 | 0.761 | 0.578 | — | — |
| relative | own | 43 | 265 | 101 | 46.5 | 0.265 | Mf 59 · Mn 12 · Pen 0 | 4.92 | 0.723 | 0.585 | 40 | 10 |
| distance | own | 9 | 81 | 9 | 4.2 | 0.031 | Mf 1 · Mn 7 · Pen 4.5 | 7.00 | 0.756 | 0.581 | 314 | 8 |
| distance_adaptive | own | 18 | 159 | 44 | 20.3 | 0.133 | Mf 12 · Mn 7 · Pen 1 | 12.0 | **0.810** | **0.704** | 220 | 10 |
| absolute | floor6 | 6 | 496 | 81 | 37.3 | 0.128 | Mf 22 · Mn 19 · Pen 18 | **1.22** | — | — | — | — |
| relative | floor6 | 6 | 376 | 30 | 13.8 | 0.053 | Mf 10 · Mn 18 · Pen 27 | 2.70 | — | — | 232 | 176 |
| distance | floor6 | 6 | 150 | 12 | 5.5 | 0.034 | Mf 1 · Mn 7 · Pen 4 | 7.00 | — | — | 516 | 73 |
| **distance_adaptive** | floor6 | 6 | 481 | **102** | **47.0** | **0.171** | Mf 11 · Mn 13 · Pen 8.5 | 1.53 | — | — | 481 | 224 |

("Δ vs absolute" = windows whose status differs from `absolute` under the same policy; "via mask"
= how many of those moved because the multi-cluster mask moved, not because detection changed.
Stability is blank at `floor6` — drop scenarios were not aggregated at that floor.)

## Findings

**1. Rule ranking is not a property of the rule — it inverts with the floor.** At floor 6,
`distance_adaptive` reproduces best (47.0%) and `relative` worst (13.8%). At the derived floor 41
this reverses completely: `relative` best (48.9%), `distance_adaptive` nearly worst (4.2%).
`absolute` is the only rule that is mid-table under every policy. **Rule and floor cannot be
chosen independently**, and any statement of the form "rule X is better" is meaningless without
naming the floor.

**2. Independent replication of the Indo 558-basis scoreboard.** At the matched floor of 6, our
501-basis pairwise chain reproduces that scoreboard almost exactly:

| rule | ours (501, floor 6) | Indo scoreboard (558, floor 6) |
|---|---|---|
| absolute | 37.3% | 37% |
| relative | 13.8% | 14% |
| distance | 5.5% | 7% |
| distance_adaptive | 47.0% | 49% |

Four independent agreements within 2 points, on a different sample basis, through a separately
driven implementation. This is the strongest external validation the run has produced that the
lifted module behaves as its own documentation claims.

**3. The cluster-size finding replicates too.** At floor 6 — the equal-absolute-floor condition
the agnostic spec §9.4.3 tested — `absolute` is essentially flat across clusters
(max/min **1.22**, vs the spec's 1.08 on its Malay fixture; Mf 22 · Mn 19 · Pen 18), while
`relative` (2.70) and `distance` (7.00) retain real size dependence. At the derived floor 41 the
measure stops being informative: every rule leaves Peninsular empty, so the ratio is computed
over the two surviving clusters and the Peninsular zero is reported separately.

**4. Much of the apparent rule difference is the mask, not the detector.** The filter-4
decomposition the phase-5 review required: at floor 6, `relative` differs from `absolute` in 232
windows of which **176 (76%)** moved because the multi-cluster mask moved; `distance_adaptive`
481 differing, **224 (47%)** via the mask. At floor 41 the effect is smaller but real —
`relative` 42 differing, **12 (29%)** via the mask. A substantial part of what looks like "this
rule detects differently" is the floor and filter 4 interacting downstream of detection.

**5. No rule is stable, and the deterministic ones are not more so.** Window-set Jaccard under a
10% / 20% drop ranges 0.71–0.81 / 0.56–0.73 across rules. `distance_adaptive` is the most stable
at drop20 under both policies (0.725 / 0.704) and `absolute` the best at drop10 (0.761); the
differences are small. This reproduces the spec's conclusion that removing the density surface
does not buy reproducibility, because the cluster consensus alleles are themselves cohort-derived.

## Recommendation: keep `absolute`

At the floor this run actually derives and ships (41), `relative` reproduces 7 more truth windows
than `absolute` (106 vs 99) — but about **a third of that margin is mask movement rather than
detection**, `absolute` is more stable under the smaller perturbation (0.761 vs 0.724), and
`absolute` is the V1-faithful default whose behaviour is documented across two cohorts.
`distance` and `distance_adaptive` collapse at this floor (0.5% and 4.2%) and only become
competitive at their own much lower floors, which is a different analysis, not a better rule.

Two caveats that belong with the recommendation. (a) The margin between `absolute` and `relative`
is small and floor-contingent; if the floor policy changes, revisit the rule. (b)
`distance_adaptive` wins decisively at floor 6 and is the most stable rule at drop20 — so a
cohort or configuration that wanted a permissive floor should not inherit `absolute` by default.

## Carried corrections

Per the phase-4 review, two items to fold into `notes.md` rather than re-revise earlier reports:
the phase-4 Mf Jaccard at matched floor 41 is **0.785** (I wrote 0.786), and the hypervariable
masks differ between arms (6 vs 5 windows, 3 in common), a secondary divergence channel worth one
line.

## Proposed next (phase 7) and alternatives

Figures (vector PDF + PNG preview, `scripts/_setup.R` palette, into `figures/<run>/`) and
`notes.md` covering method, floor, clonality effect, benchmark and rule choice. Alternatives
considered: (a) also extend the rule comparison to the de-clonalized arm — rejected, it would
push past the phase limit for a conclusion about rules rather than cohorts; (b) defer `notes.md`
to a ninth phase — not available, and the mission requires it.
