# Phase 4 — aggregation at the derived floors, both arms

**Scope (reviewer's phase-3 spec).** Run `introgression_aggregate.R` per arm at floor 41 (full) /
40 (de-clonalized) with config-sourced gene-family filters; produce the filtered window table,
per-pair and per-sample tables, and the filter-audit chain; call out Mf↔Mn; quantify the
clonality effect **on the filtered sets**; report Peninsular's zero rows as-is.

**Done-condition.** Both arms' filtered window + per-sample tables non-empty; Mf↔Mn window and
sample counts per direction; clonality effect via window-set Jaccard per cluster and per-sample
call-count deltas over shared samples; every number traceable. — **Met.**

---

## What was run

Rerun cleanly into `aggregate/` (the reviewer's preferred option) rather than copying from the
sweep work dirs — 1.9 s and 1.7 s respectively, so there was no reason to reuse. The outputs
**reproduce the sweep's rows exactly** (full N=41: 24,658 calls / 263 windows / 471 samples;
declonal N=40: 22,946 / 254 / 456), which is an independent confirmation that the sweep and the
aggregation agree.

```
scripts/introgression_module/introgression_aggregate.R
  --clusters <arm table>  --metadata inputs/samples.tsv  --contig-map inputs/contig_map.tsv
  --fai data/gadi/fasta/strain_A1_H.1.Icor.fasta.fai
  --gff data/gadi/gff/strain_A1_H.1.Icor.gff3  --gene-family-filters "SICA,KIR"
  --window-size 10000 --min-samples-per-window 2
  --per-cluster-min-pct 0 --per-cluster-min-samples <41|40>
  --out-dir results/<run>/aggregate/<arm>   calls/<arm>/*.tsv
```

**GFF mapping resolved cleanly:** 14 of 16 contigs mapped, all by exact **name** — no length or
nearest-length fallback was needed, because this project's GFF and reference share the
`ordered_PKNH_NN_v2` naming. The 2 unmapped are the mitochondrial and apicoplast contigs,
correctly excluded from the nuclear set. This is the V1 defect the module fixed (agnostic spec
§9.1.5), and on our inputs it simply does not arise. 1,298 SICA/KIR features across all 14
contigs → a 384-window gene-family mask.

## Outputs

`results/<run>/aggregate/{full,declonal}/` — `introgressed_windows_filtered.tsv`,
`windows_by_cluster.tsv`, `windows_across_chrom.tsv`, `intro_per_sample_summary.tsv`,
`window_sample_counts_raw.tsv`, `average_windows_for_clusters.tsv`, `filter_audit.tsv`,
`gene_family_masked_windows.tsv`, `hypervariable_masked_windows.tsv`, `gff_contig_map.tsv`,
`aggregate.log`. **All of this lives only under gitignored `results/` — the per-sample tables
carry sample IDs and are never tracked.**

## Key numbers

**Filter-audit chain** (`aggregate/<arm>/filter_audit.tsv`):

| step | full: calls | windows | samples | declonal: calls | windows | samples |
|---|---|---|---|---|---|---|
| raw (all pairs) | 53,615 | 1,190 | 501 | 50,970 | 1,205 | 484 |
| 1. dataset n > 2 | 53,584 | 1,173 | 501 | 50,930 | 1,187 | 484 |
| 2. per-cluster floor | 26,275 | 279 | 471 | 24,686 | 269 | 456 |
| 3. gene-family mask (384 win) | 25,520 | 269 | 471 | 23,750 | 259 | 456 |
| 4. hypervariable (6 / 5 win) | **24,658** | **263** | **471** | **22,946** | **254** | **456** |

The support floor does essentially all the work: it removes 51% of calls and 76% of windows. The
dataset-level `n > 2` filter removes 31 calls out of 53,615 — the same "doing almost nothing"
result the spec records for the Indo cohort (§9.2).

**Per cluster at the floor** (`windows_by_cluster.tsv`):

| arm | Mf windows | Mf samples | Mn windows | Mn samples | Peninsular |
|---|---|---|---|---|---|
| full (N=41) | 249 | 366 | 14 | 105 | **0 windows** |
| declonal (N=40) | 244 | 359 | 10 | 97 | **0 windows** |

**Headline pair Mf↔Mn, per direction** (`introgressed_windows_filtered.tsv`, `PAIR == "Mf__Mn"`):

| arm | direction | calls | windows | samples |
|---|---|---|---|---|
| full | Mf_like_Mn | 10,719 | 222 | 366 |
| full | Mn_like_Mf | 236 | 11 | 101 |
| declonal | Mf_like_Mn | 10,266 | 216 | 359 |
| declonal | Mn_like_Mf | 214 | 9 | 93 |

The pair is **strongly asymmetric**: Mf-carrying-Mn-signal outnumbers the reverse ~45:1 in calls
and ~20:1 in windows. Any "Mf↔Mn intermixing" wording should say which direction carries it.

**Chromosomal distribution** (`windows_across_chrom.tsv`): the Mf signal concentrates on
chr11–14 — 36 + 43 + 30 + 35 = **144 of 249 windows (58%)** on four of fourteen chromosomes.
This matches the legacy script's margin note that chr12/13/14 carry most of the genome-wide
signal (agnostic spec §9.2), reached here independently.

**Samples retained**: every cluster member that survives filtering is *all* of them — Mf 366/366,
Mn 105/105 (full). Extending the reviewer's phase-2 point: saturation is not just a raw-detection
property, it survives the entire filter chain.

## The clonality effect, measured on the filtered sets — and a correction to phase 2

Comparing full@41 with declonal@40 conflates two things, so I decomposed it
(`aggregate/` + the sweep work dirs `_sweep_full/N40`, `_sweep_declonal/N41`):

| comparison | what it isolates | Mf Jaccard | Mn Jaccard |
|---|---|---|---|
| full@40 vs full@41 | **floor alone** | 0.969 | 0.933 |
| full@40 vs declonal@40 | **clonality alone** | 0.776 | 0.667 |
| full@41 vs declonal@41 | **clonality alone** | 0.786 | 0.500 |
| full@41 vs declonal@40 | as reported (both) | 0.780 | 0.714 |

**The floor difference is negligible; de-clonalization is what moves the filtered window sets.**
Per-sample filtered call counts over the 456 shared samples: mean −2.61, median −2, range −12 to
+10; 340 decreased, 37 unchanged, **79 increased**.

**This corrects my phase-2 framing.** I wrote that de-clonalization is "a clean subtraction at
the sample level". That is true for sample *status* — who is called at all — but **not** at the
call level. Among retained samples in the `Mf__Mn` pair, raw (sample, window) calls agree at
Jaccard **0.883**: 20,303 shared, 1,317 only-full, 1,369 only-declonal — and **454 of 456**
retained samples have at least one changed call. Removing 17 redundant genotypes re-derives the
cluster consensus and re-shapes the contours, so retained samples keep being called but partly in
*different windows*. That is why the filtered window sets diverge (0.78) far more than the raw
window sets did (0.99): the filter chain aggregates per-window support, and ~12% of calls moving
is enough to reshuffle which windows clear a floor of 41.

I also checked the obvious simpler explanation — windows sitting just above the floor being
tipped below it — and it is **not** sufficient. Of the 37 Mf windows present at full@40 but absent
at declonal@40, only 15 had support ≤ 50; 22 had support ≥ 51, and the median lost window had 59
of 366 samples. Only 7 Mf samples were removed, so a drop from 59 to below 40 cannot come from
the removal itself. The consensus/contour recomputation is doing it, as the call-level numbers
above show directly.

## Findings

1. **Peninsular contributes zero windows to the filtered result in both arms**, as derived in
   phase 3. Reported as-is; not tuned around. The aggregate log reads "2 clusters" for this
   reason. The NEEDS-JACOB stands unchanged.
2. **The headline survives and is quantified, but it is one-directional and Mn-thin** — 222 vs 11
   windows in the full arm. The Mf side is robust; the Mn side is 11 windows and 9 after
   de-clonalization.
3. **The filtered window set is the cohort-sensitive layer.** Sample status is stable, raw window
   sets are stable (0.99), but filtered window sets are not (0.78) — the floor amplifies small
   call-level movements. This is the single most important caveat for the benchmark in phase 5:
   we are comparing a cohort-sensitive object against a published cohort-sensitive object.
4. **The gene-family filter behaved as designed** and the V1 contig-mapping defect does not arise
   on this project's inputs (exact-name mapping, all 14).

## Proposed next (phase 5) and alternatives

Coordinate benchmark vs the poster's 217 windows — CHROM + window start, never WINDOW id — for
both arms, reporting % reproduced alongside the Indo 558-basis 64% ceiling and the 501 ∩ 558
overlap, per the phase-0 mitigation. Alternatives considered: (a) benchmark only the full arm as
the closer match to the poster's cohort — rejected, the mission requires both arms; (b) benchmark
the unfiltered window sets to get a higher overlap — rejected, that compares different objects,
since the truth is a filtered result.
