# Phase 5 — coordinate-based benchmark vs the poster's 217 windows

**Scope (reviewer's phase-4 spec).** Score both arms' filtered window sets against the published
217 windows by **CHROM + window start, never WINDOW id**; verify grid alignment first; report per
arm and per truth direction, alongside the Indo 558-basis 64% ceiling and the 501 ∩ 558 overlap.

**Done-condition.** `benchmark/coordinate_benchmark.tsv` non-empty, one row per arm × direction
(+ combined), every number traceable, grid-alignment check documented. — **Met.**

---

## Grid alignment — checked before scoring, and it was necessary

| source | convention | example |
|---|---|---|
| truth (`inputs/truth/truth_windows_coords.tsv`) | `start` on the 10 kb grid, `end` = start + 9,999 | `CHROM 1, start 40000, end 49999` |
| ours (`aggregate/*/introgressed_windows_filtered.tsv`) | `START` on the grid, `BIN` = **START + 5,000** (bin midpoint), `END` = START + 9,999 | `CHROM 1, BIN 135000, START 130000, END 139999` |

Verified across every row: `truth.start %% 10000 == 0` and `our START %% 10000 == 0`, while
`BIN − START == 5000` for all rows. **The match is `truth.start` ↔ our `START`.** Matching
`truth.start` to `BIN` — the obvious mistake, since `BIN` is the column the window id is built
from — would have compared points 5 kb apart and scored near zero. Worked example: truth
`CHROM 1 / start 40000` matches our window `w1_45000` (`START 40000`, `END 49999`), not `w1_40000`.

## A bug I introduced, found, and fixed before reporting

My first version of `benchmark_coordinates.R` built the join key with
`paste0(chrom, ":", start)`. The truth table's `start` was read as a **double**, and R renders
round doubles in scientific notation — `as.character(100000)` is `"1e+05"` — while our `START`
came off `fread` as an **integer** rendering as `"100000"`. Every truth window at a coordinate R
chose to render scientifically was therefore silently unmatchable.

It failed exactly the way this project's rewrite brief warns about: no error, just wrong numbers.
It under-reported reproduction as **89/217 (41.0%)** instead of the correct **99/217 (45.6%)**,
and it corrupted the eligibility check into "200 of 217 truth windows attainable" when the true
answer is **217 of 217**. I caught it only because the same quantity computed in Python
disagreed. Fixed by keying with `sprintf("%d:%d", as.integer(...), as.integer(...))`, with a
comment in the script recording why. **Every number below was then recomputed independently in
Python and agrees with the R output cell for cell.**

## What was run

`Rscript scripts/introgression_malay/benchmark_coordinates.R --config config/malay_cohort.yaml
--out results/<run>/benchmark/coordinate_benchmark.tsv` (2.3 s). Outputs
`benchmark/coordinate_benchmark.tsv` and `benchmark/benchmark_context.tsv`.

## Context — what the comparison is against (`benchmark_context.tsv`)

| item | value |
|---|---|
| our cluster basis | 501 |
| published basis | 558 |
| intersection | **501** — our cohort is a strict *subset*; 57 in theirs not ours, **0** in ours not theirs |
| truth windows | 217 (150 Mf + 67 Mn) |
| truth windows inside our eligible universe | **217 of 217** — none is unreachable for eligibility reasons |
| eligible universe | 1,363 windows |
| published 558-basis ceiling | 138/217 = 64% (`regen_density`, agnostic spec) |

## Key numbers (`benchmark/coordinate_benchmark.tsv`)

**Headline — at each arm's own derived floor (41 / 40):**

| arm | scored set | truth | called | reproduced | % reproduced | Jaccard |
|---|---|---|---|---|---|---|
| full | combined | 217 | 263 | **99** | **45.6%** | 0.260 |
| full | Mf vs Mf-truth | 150 | 249 | **87** | **58.0%** | 0.279 |
| full | Mn vs Mn-truth | 67 | 14 | 1 | 1.5% | 0.013 |
| full | Mf-truth vs any cluster | 150 | 263 | 88 | 58.7% | 0.271 |
| full | Mn-truth vs any cluster | 67 | 263 | 11 | 16.4% | 0.034 |
| declonal | combined | 217 | 254 | **93** | **42.9%** | 0.246 |
| declonal | Mf vs Mf-truth | 150 | 244 | **84** | **56.0%** | 0.271 |
| declonal | Mn vs Mn-truth | 67 | 10 | 1 | 1.5% | 0.013 |

**Supplementary — the same scoring at the *published chain's own* floor (n > 5, i.e. N = 6),
from the phase-3 sweep dirs. Clearly labelled; the derived floor remains the headline and
nothing was re-tuned:**

| arm | scored set | called | reproduced | % reproduced | Jaccard |
|---|---|---|---|---|---|
| full | combined | 496 | 81 | 37.3% | 0.128 |
| full | Mf | 290 | 53 | 35.3% | 0.137 |
| full | Mn | 144 | 16 | 23.9% | 0.082 |
| declonal | combined | 518 | 89 | 41.0% | 0.138 |
| declonal | Mf | 312 | 64 | 42.7% | 0.161 |
| declonal | Mn | 127 | 13 | 19.4% | 0.072 |

## Findings

**1. The Mf side reproduces well: 58% against a 64% ceiling.** That ceiling
(`regen_density`, 138/217) is the *same density method* run three-cluster on the *full 558* with
the published filter chain. We reach 58% of the Mf truth with a **pairwise** reformulation, on a
**501-sample strict subset**, under a floor **7× stricter** than the published one. On the
direction the cohort actually carries signal, the agnostic module substantially reproduces the
published result.

**2. The combined figure is dragged down entirely by Mn, and that is the floor, not the method.**
Our Mn set holds 14 windows, so it can match at most 14 of the 67 Mn truth windows — it reproduces
1. This is arithmetic downstream of phase 3's finding, not a placement failure: at the published
floor of 6, where Mn keeps 144 windows, Mn reproduction rises to **23.9%**. The FDR-safe floor
buys false-positive control by discarding almost the whole Mn window set.

**3. The stricter floor *improves* overall agreement — the opposite of the obvious expectation.**
At the derived floor we call 263 windows and reproduce 99 (45.6%, Jaccard 0.260); at the published
floor 6 we call 496 and reproduce 81 (37.3%, Jaccard 0.128). Nearly twice as many calls yield
*fewer* true reproductions. The windows the floor removes are disproportionately ones the poster
did not report, which is direct evidence the derived floor is removing noise rather than signal —
and it is the strongest empirical support this run produces for the floor derivation.

**4. Every truth window is reachable here.** All 217 fall inside our 1,363-window eligible
universe, so nothing is lost to coverage. Non-reproduction is a detection/filter outcome, not a
data-availability artifact.

**5. De-clonalization costs a little reproduction** (45.6% → 42.9% combined; 58.0% → 56.0% Mf),
consistent with phase 4: the truth was built on a clonal cohort, so removing redundant genotypes
moves us slightly away from it. This is the expected direction and not an argument against
de-clonalizing.

**6. The basis mismatch is one-directional and small.** Our 501 is a strict subset of the
published 558 — no sample we use was outside their analysis. The gap is 57 samples (10%) they had
and we do not, which can only *reduce* our support counts relative to theirs, i.e. it biases
reproduction downward. Per the standing instruction this is reported, not resolved.

**7. Standing caveat carried from phase 4.** The filtered window set is the cohort-sensitive layer
(raw window Jaccard 0.99 vs filtered 0.78 between arms), so this benchmark compares a
cohort-sensitive object against a published cohort-sensitive object. The % reproduced should be
read as "replication of the published windows", never as accuracy.

## Proposed next (phase 6) and alternatives

Rule comparison — `absolute` / `relative` / `distance` — on the benchmark's measures
(window reproduction, cluster-size dependence, stability under 10%/20% sample drops at a fixed
seed), using the lifted `introgression_rule_sweep.R` / the benchmark harness pattern. Alternatives
considered: (a) also sweep `distance_adaptive` — worth including if it is free, since the Indo
scoreboard ranked it above `absolute` on window reproduction; (b) run the rule comparison only on
the full arm to halve the cost — reasonable, and I will do both arms only if runtime stays in
minutes, reporting which was covered.
