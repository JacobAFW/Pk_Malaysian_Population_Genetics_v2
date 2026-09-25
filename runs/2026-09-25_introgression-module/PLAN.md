# PLAN — agnostic introgression module on the Large-Malay cohort

Run: `runs/2026-09-25_introgression-module/` · branch `auto/2026-09-25_introgression-module`
Mission: `MISSION.md` · protocol: `AUTONOMY_PROTOCOL.md` · phase limit 8 · 3 threads.

## Facts established before planning (all verified on disk)

- All 8 mission inputs present. Genotype `data/processed/introgression/hmmIBD.tsv`
  = 61,951 variants × 755 sample columns. Labels `data/processed/cluster_maf/cluster_labels.tsv`
  = **501 samples: Mf 366 · Mn 105 · Peninsular 30**; all 501 are columns of the genotype table.
- Poster truth for the benchmark is held locally: `data/raw/Introgression/introgressed_windows_filtered.tsv`
  (the 3,942-row / 217-window table). Reference `.fai` exists (`data/gadi/fasta/strain_A1_H.1.Icor.fasta.fai`),
  so no bed→fai reconstruction is needed.
- Env `vvg-box-pixi` has R 4.5.3 with data.table/dplyr/tidyverse/sp/ggplot2. The module's R
  scripts use a hand-rolled `--flag value` parser — **`optparse` is missing but not needed**.
- Indo module (read-only) already ships a Malay benchmark harness
  (`scripts/R/introgression_benchmark_malay.R`, 4 rules × 3 scenarios × 3 pairs in ~4.7 min).
  Its fixture is a 2026-08-17 copy of *this* cohort, so the method is known to run at this scale.

## Decisions taken here (flagged for the reviewer)

1. **Driver = direct `Rscript`, not Snakemake.** The module's scripts have complete CLIs and the
   Indo benchmark already drives them directly; a local Snakemake target would add a dependency
   (snakemake in vvg-box is unverified) for 3 pairs. Alternative considered: lift
   `workflow/rules/05_introgression.smk` — rejected as cost without benefit at K=3.
2. **Cluster basis = `cluster_labels.tsv` (501)**, as the mission names it. This is *not* the
   558-sample set the poster truth was built on — a known open item (`DECISIONS.md`: ADMIXTURE
   Q basis). Recorded as a benchmark caveat, not resolved here.
3. **Rule comparison via a lifted copy of the benchmark harness** (`introgression_benchmark_malay.R`
   → local `introgression_benchmark_cohort.R`), input-prep swapped to our labels + `.fai`. It
   already computes the distance table once per pair and calls every rule against it, plus the
   drop-10%/20% stability and the cluster-size test. Alternative: reimplement — rejected (§ mission
   "never edit in place"; a lifted copy satisfies both).
4. **Focal-subgroup test is NOT run** (mission §5, resolved N/A).

## Phases

| # | Scope | Done when | Outputs |
|---|---|---|---|
| 1 | Lift module into `scripts/introgression_module/`; build the 5 derived input tables (clusters, metadata, contig_map, fai, de-clonalized keep-list); write `config/malay_cohort.yaml`; smoke-test one pair (Mf__Mn) end-to-end | one pair's call table exists with non-zero rows; every input's provenance recorded; clonal collapse reports n kept/dropped per cluster | `scripts/introgression_module/**`, `config/malay_cohort.yaml`, `results/2026-09-25_introgression-module/inputs/*.tsv`, phase report |
| 2 | Detection, `absolute` @5e-4, all 3 pairs × both arms (full 501, de-clonalized) | 6 per-pair call tables written and non-empty | `results/<run>/calls/{full,declonal}/*.tsv` |
| 3 | Floor derivation on THIS cohort: permutation null (1,000 reps, stated seed) + explicit FDR target, per cluster, both arms; stability sweep around the chosen N | floor table gives per-cluster null mean/p95/p99/max, FDR-safe N, expected-null-windows and observed windows; sweep table spans N | `results/<run>/floor/floor_derivation_{full,declonal}.tsv`, `floor_stability_sweep_*.tsv` |
| 4 | Aggregate under the derived floor: 4-filter chain, both arms → per-pair window tables + per-sample tables; headline pair Mf↔Mn called out | `introgressed_windows_filtered.tsv` + per-sample counts exist for both arms; clonality effect quantified | `results/<run>/aggregate/{full,declonal}/*.tsv` |
| 5 | **Coordinate** benchmark vs the poster's 217 windows (CHROM + window start, never WINDOW id), both arms → % reproduced | benchmark table gives shared/217, Jaccard, and %-reproduced per arm, traceable to the truth file | `results/<run>/benchmark/coordinate_benchmark.tsv` |
| 6 | Rule comparison `absolute`/`relative`/`distance` on the benchmark measures: window reproduction, cluster-size dependence (rho + max/min), stability under 10%/20% sample drop (fixed seed) | scoreboard covers 3 rules × the measures; a recommended rule is stated with reasoning | `results/<run>/benchmark/rule_scoreboard.tsv` |
| 7 | Figures (vector PDF + PNG preview, `scripts/_setup.R` palette) + `notes.md` | figures render in the new versioned folder; notes cover method, floor, clonality effect, benchmark, rule choice | `figures/2026-09-25_introgression-module/*.{pdf,png}`, `results/<run>/notes.md` |
| 8 | Reserve — revision of whichever phase the reviewer returns REVISE on, or run close-out | reviewer PROCEED, or stop condition met | `RUN_SUMMARY.md` |

## Standing constraints applied

- Indo `agnostic/` is read-only: copied out, never written to. `data/raw/`, `data/gadi/`,
  existing `figures/`, and `main` are untouched. All new outputs are in new versioned folders.
- `results/` and `figures/` are gitignored; `runs/` is tracked — **reports carry counts only,
  no sample IDs**. commit-guard runs before every commit.
- No installs; a missing package is a blocked action. Local only, no pushes.
- Each phase's numbers cite the file they can be checked against.

## Known risks

- The 501-vs-558 cohort mismatch will depress the coordinate benchmark for reasons unrelated to
  the method. Will be reported as a caveat with the overlap actually measured, not tuned around.
- Peninsular n=30 (n smaller after de-clonalization) — an FDR-safe global floor set by Mf may
  again leave Peninsular with zero windows (the Indo §9.5 finding). If so it is reported as a
  finding and tagged NEEDS-JACOB, not tuned around.
