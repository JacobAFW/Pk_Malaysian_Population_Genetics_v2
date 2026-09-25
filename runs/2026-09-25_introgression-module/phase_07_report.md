# Phase 7 — figures and notes.md

**Scope (reviewer's phase-6 spec).** Publication figures (vector PDF + PNG preview) into a new
versioned `figures/<run>/` using the `scripts/_setup.R` palette, reading only existing TSVs; plus
`results/<run>/notes.md` covering method, floor, clonality, benchmark and rule choice, carrying
the corrections register and all NEEDS-JACOB items; and close the `agg_floor6` reproducibility gap.

**Done-condition.** Valid PDF+PNG pairs; notes.md complete with every headline number citing its
file; `score_rules.sh` regenerates `agg_floor6`; no sample IDs in figures or notes. — **Met.**

---

## What was run

1. `scripts/introgression_malay/score_rules.sh` — added the `floor6` branch (the reviewer's
   finding 4: the floor-6 aggregations underpinning the run's strongest validation existed on
   disk but no committed script produced them). **Verified by regeneration**: moved the existing
   output aside, re-ran, and the new `introgressed_windows_filtered.tsv` is byte-identical
   (md5 `30a1572d…`). The headline validation is now reproducible from the committed code.
2. `scripts/introgression_malay/make_figures.R` — 5 figures, no analysis, reads only existing
   TSVs under `results/<run>/`.
3. `results/2026-09-25_introgression-module/notes.md`.

## Outputs

`figures/2026-09-25_introgression-module/` — 5 × (PDF + PNG), all confirmed valid
(`file`: "PDF document, version 1.4, 1 pages"; PNGs 2700 px wide at 300 dpi). Existing `figures/`
untouched; this is a new folder.

| figure | content |
|---|---|
| `fig1_floor_derivation` | (a) windows surviving each floor per cluster, both arms, floor marked; (b) null support per cluster vs the derived floor |
| `fig2_floor_filter4_interaction` | multi-cluster mask size and surviving windows vs N, both floors marked |
| `fig3_headline_mf_mn` | (a) Mf↔Mn directional asymmetry; (b) per-chromosome distribution, Peninsular explicitly at zero |
| `fig4_clonality_effect` | (a) raw vs filtered Jaccard per cluster; (b) per-sample window-count change |
| `fig5_benchmark_and_rules` | (a) reproduction vs the 64% ceiling; (b) rule ranking by floor policy with the published scoreboard marked |

## Figure defects I found by looking at the rendered output

Each was caught by inspecting the PNGs rather than trusting the code, and each is fixed:

1. **Fig 4a compared unlike things.** The raw-layer Jaccards were per **pair** (Mf__Mn 0.993,
   Mn__Peninsular 0.930) while the filtered-layer ones were per **cluster**, both plotted under a
   single "Cluster" legend — a mislabelling that would have been read as a per-cluster
   comparison. Recomputed the raw layer per cluster from the call tables: **Mf 0.969, Mn 0.966,
   Peninsular 0.844** against filtered Mf 0.776 / Mn 0.667. The qualitative conclusion (the
   filtered layer is the sensitive one) survives and is now measured consistently.
2. **Fig 5b's comparison markers sat over the wrong bar** — plotted at the rule's centre rather
   than on the floor-6 bar they annotate. Now positioned on the matching dodged bar, which is
   also what makes the agreement visible.
3. **Facet order was alphabetical**, putting the de-clonalized arm left of the full arm. The full
   arm is the reference and now comes first.
4. **Em-dashes were substituted in the PDF device** (`_setup.R` uses base `pdf()` because cairo
   needs X11). All figure text is now ASCII, so the vector output matches the preview.
5. **Fig 4b's subtitle over-claimed its sample set** — the 456 samples are Mf/Mn only, because
   Peninsular retains no filtered window. Stated explicitly rather than leaving the absence
   unexplained.

## notes.md

Ten sections: method and lift provenance; the derived floor with the Peninsular incompatibility;
clonality (including the named-input substitution and the status-vs-call-level distinction);
results at the floor with the Mf↔Mn asymmetry; the coordinate benchmark with the 64% ceiling,
the 501⊂558 basis caveat and the floor×filter-4 mechanism; rule choice with the ranking-inversion
finding and the `absolute` recommendation; a **corrections register** of all seven superseded
numbers/claims; scope and coverage; the four NEEDS-JACOB items; and the figure index.

Every headline number cites its file. The register exists so no superseded figure reaches a
draft — including the two the reviewer required (phase-4 Jaccard **0.785** not 0.786; the phase-6
stability superlative, where `distance_adaptive` leads at both drop levels and `absolute` leads
only `relative`).

## Checks

- No sample identifiers in any figure (per-sample panels are histograms, never labelled points)
  or in `notes.md` — scanned file-by-file, zero hits.
- `notes.md` lives under gitignored `results/`; figures under gitignored `figures/`. Only the
  scripts and this report are tracked.
- Nothing recomputed: no detection, floor derivation, aggregation or benchmark was re-run. The
  one write outside `figures/<run>/` was the `score_rules.sh` floor6 fix and its verification.

## Blocked action (protocol §5)

The floor6 regeneration check left a scratch copy at
`results/<run>/rules/full/absolute/agg_floor6_prev/`. The guard denies deletes outside scratch
dirs, so it remains. It unlocks nothing and is inside gitignored `results/`; recorded rather than
worked around. No `NEEDS-JACOB`.

## Proposed next

Close the run per protocol §6: `RUN_SUMMARY.md` with the NEEDS-JACOB queue first, then what
changed, the low-confidence reviewer decisions, and a suggested next run (de-clonalized rule
sweep; the Peninsular own-scale null; the introgression × environment association). This is
phase 7 of 8, so the summary fits inside the limit.
