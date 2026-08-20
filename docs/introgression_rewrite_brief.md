# Introgression pipeline — known pitfalls brief (for the rewrite)

**Purpose.** The old windowed-introgression script (`find_introgression_updated_framework.R`)
had a silent, version-fragile bug *and* a deeper method limitation. This brief captures both
so the new data-agnostic pipeline doesn't re-inherit them. Full diagnosis:
`results/stage_b/19_diagnosis_report.md`.

## What went wrong

### 1. The bug — positional column indexing on ggplot2 contour internals
The old script assigned samples to introgressed regions by (a) building a 2-D density surface
with `ggplot2::geom_density_2d` over the sample-distance point cloud, (b) pulling contour-polygon
coordinates **by positional index** out of ggplot's internal layer data (`Mf_contours[3]` for x,
`Mf_contours[4]` for y), and (c) classifying points with `sp::point.in.polygon`.

ggplot2 4.0.3 **reordered** the contour-layer columns, so `[3]` is now `y` and `[4]` is `piece`
(an integer segment ID). The polygon coordinates were therefore garbage and `point.in.polygon`
returned effectively random in/out calls. Symptom: ~**6.3× too many** introgressed windows
(24,980 vs 3,942 truth), varying by ggplot2 version. **It failed silently — no error, just wrong
numbers**, and the pre-fix output even had *coincidental* overlap with the truth windows, which
masked it. (Immediate fix in the ported script: use named columns `c$x`, `c$y` — but see rule 2.)

### 2. The method limitation — density contours are cohort-sensitive by construction
Even with the bug fixed, the contour approach **cannot reproduce across cohorts**:
`geom_density_2d` builds the surface from the exact sample point cloud, so any change in cohort
composition or coordinates shifts the `level > 5e-4` polygons → different windows. In replication
this left an irreducible ~2.25× window-count gap and only ~19% window-ID overlap **despite the
sample-level result being exact**. The "which samples show introgression" answer is robust; the
"which genomic windows" answer is an artifact of the density-surface method.

## Design rules for the new version

1. **Deterministic window assignment.** Replace `geom_density_2d` + `point.in.polygon` with a rule
   that depends only on the data, not a fitted density surface — e.g. per-window allele-match
   fraction between a sample and each cluster, thresholded. Same input → same windows, across
   cohorts and runs.
2. **Never positional-index a library's internal data structures** (ggplot layer data, contour
   output, etc.). Use named columns and assert them; better, don't build analysis logic on
   plotting internals at all. *This was the actual root-cause bug.*
3. **FDR / multiple-testing correction** (BH) on the window-significance regressions. The old
   pipeline used raw p<0.05 across ~150 windows × 4 environmental terms × interactions →
   over-detection (the 34-vs-24-windows finding). Report both raw and adjusted.
4. **Version-robust dplyr.** `reframe()` for multi-row-per-group output (not `summarise()`);
   `consecutive_id()`/`dense_rank()` for group IDs (not `group_indices()`). Avoid deprecated idioms
   whose behaviour changes across versions.
5. **Join by ID, never by row position.** The coord file, keep-list and `.fam` were in three
   different orders in the old data. (Also: no sample IDs hard-coded in the script — hold
   exclusion/drop lists in gitignored external files; the repo is public.)
6. **Add a synthetic-cohort ground-truth test** — a small simulated dataset with known introgressed
   windows — so correctness is checkable deterministically and diffable across runs (old outputs
   were cohort-specific and couldn't be compared cleanly).

## Verification checklist for the new script
- [ ] No positional indexing of ggplot/contour/library-internal objects.
- [ ] Window assignment is deterministic (run twice on same input → identical window set).
- [ ] p-values BH-adjusted; raw + adjusted both reported.
- [ ] All sample/coord joins are by ID.
- [ ] No sample IDs inline (externalised to gitignored files).
- [ ] Passes the synthetic ground-truth case.
