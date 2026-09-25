# Phase 1 — lift the module, assemble the inputs, smoke-test one pair

**Scope (reviewer's phase-0 spec).** Lift the introgression module into this project with
provenance; build the derived input tables (clusters, metadata, contig map, fai, de-clonalized
drop-list); copy the missing `mn_windows.tsv` coordinate truth and verify it covers all 217
filtered WINDOW ids; write the cohort config; smoke-test pair `Mf__Mn` end-to-end.

**Done-condition.** `Mf__Mn` call table exists with non-zero rows; input provenance complete;
clonal collapse reported as counts per cluster; nothing written under the read-only module path,
`data/raw/`, or `data/gadi/`. — **Met.**

---

## What was run

| step | command |
|---|---|
| lift | `cp -n` of 8 `introgression_*.R` + `LICENSE` from the read-only module `scripts/R/` → `scripts/introgression_module/`; md5 of every copy recorded against its source |
| prep | `Rscript scripts/introgression_malay/prepare_inputs.R --config config/malay_cohort.yaml --out-dir results/<run>/inputs` |
| smoke | `Rscript scripts/introgression_module/introgression_pair.R --genotype-table data/processed/introgression/hmmIBD.tsv --clusters <run>/inputs/admix_clusters.tsv --pair Mf__Mn --window-size 10000 --min-snps 5 --detection-rule absolute --contour-level-other 5e-4 --contour-level-own 5e-4 --out <run>/calls/full/Mf__Mn.tsv` |

Source module commit `d6a788e` — all 9 lifted files md5-identical to source
(`scripts/introgression_module/PROVENANCE.md`). The module tree was not written to.

## Outputs

| path | contents |
|---|---|
| `scripts/introgression_module/` | 8 lifted module scripts + `LICENSE` + `PROVENANCE.md` |
| `scripts/introgression_malay/prepare_inputs.R` | the cohort adapter (this project's only new analysis code so far) |
| `config/malay_cohort.yaml` | every cohort-specific path and parameter |
| `results/<run>/inputs/` | `admix_clusters.tsv`, `samples.tsv`, `contig_map.tsv`, `sample_missingness.tsv`, `declonalization_audit.tsv`, `declonalization_counts.tsv`, `exclude_{full,unique}.txt`, `{all,unique}_genotypes.txt`, `input_provenance.tsv` |
| `results/<run>/inputs/truth/` | `mf_windows.tsv`, `mn_windows.tsv`, `introgressed_windows_filtered.tsv`, `truth_windows_coords.tsv` |
| `results/<run>/calls/full/Mf__Mn.tsv` | the smoke-test call table |

`results/` and `figures/` are gitignored; only `scripts/`, `config/` and `runs/` are tracked.

## Key numbers

| number | value | check against |
|---|---|---|
| clustered samples | 501 — Mf 366 · Mn 105 · Peninsular 30 | `inputs/admix_clusters.tsv` |
| genotype table | 61,951 variants × 755 sample columns; all 501 present | `input_provenance.tsv` (md5 of `hmmIBD.tsv`) |
| contigs | 14 nuclear (2 non-nuclear dropped); CHROM codes match the genotype table exactly | `inputs/contig_map.tsv` |
| per-sample missingness | median F_MISS 0.0008, range 0.0000–0.2173 over 61,951 sites | `inputs/sample_missingness.tsv` |
| clonal pairs ≥ 0.98 | 20 in the full IBD table; **19** with both members clustered | `inputs/declonalization_audit.tsv` |
| clonal groups | **15** (13 of size 2, 2 of size 3) covering 32 samples; 0 straddle a cluster | `inputs/declonalization_audit.tsv` |
| de-clonalized arm | 501 → **484** (17 dropped): Mf 366→359 · Mn 105→97 · Peninsular 30→28 | `inputs/declonalization_counts.tsv` |
| coordinate truth | 150 Mf + 67 Mn = 217 windows, disjoint, covering **217 of 217** filtered WINDOW ids | `inputs/truth/truth_windows_coords.tsv` |
| smoke test `Mf__Mn` | 473,977 eligible sample-windows → **22,371 calls / 1,040 windows / 471 samples** (Mf_like_Mn 16,450 · Mn_like_Mf 5,921) | `calls/full/Mf__Mn.tsv` (22,371 data rows) |
| smoke-test cost | 41 s wall, 5.3 GB peak RSS | `/usr/bin/time -l` |

At ~41 s and 5.3 GB per pair, the six detection runs in phase 2 (3 pairs × 2 arms) cost ~4 min —
far inside the 2 h/step budget. Threads capped at 3 (`setDTthreads(3L)`).

## What went wrong / had to change

**1. The named clonal inputs are unusable as given — and would have failed silently.**
`MISSION.md` §3 names `data/raw/hmmIBD/IBD_Clonal_transmission_miss5.tsv` and
`Mn_clonal_transmission.tsv` as the clonal sets. They are in the **field-collection id
namespace**; the genotype table and cluster labels use sequencing ids. The intersection is
**0 of 21** named samples (`prepare_inputs.R` reports this explicitly). The first run of the
prep script accordingly produced a de-clonalized arm identical to the full arm — 0 dropped —
which is exactly the kind of silent no-op the introgression rewrite brief warns about.

**Resolution, without inventing anything.** The de-clonalized arm is built from this project's
own full pairwise hmmIBD output, `data/gadi/validation/hmmIBD/Pk.hmm_fract.txt` (155,403 pairs
over 558 samples, all in the genotype namespace, covering all 501 clustered samples). The
threshold is not a new choice: `scripts/14b_clonal_networks.R` states this project's rule as
`fract_sites_IBD >= 0.98`, and the named tables' own minimum value (0.98672) is consistent with
it. The cut is well separated — 20 pairs at ≥ 0.98, only **3** in [0.90, 0.98) — so it is not
perched on a shoulder. Both the source and the threshold are declared in
`config/malay_cohort.yaml`; the named tables are retained as a cross-check and reported on every
run. Alternative considered: declare the de-clonalized arm blocked and NEEDS-JACOB — rejected,
because the project holds a usable source and the threshold is the project's own, so blocking
would drop half the mission's deliverable for no evidentiary gain.

**2. No `.smiss`/`.imiss` exists anywhere in the project.** The module's representative rule
(`docs/clonality.md`) is "lowest per-sample missingness, alphabetical tie-break", degrading to
alphabetical-only without a `.smiss`. Rather than degrade, F_MISS is computed directly from the
genotype table (missing = `-1`) — the same quantity PLINK reports, from the same calls,
deterministic and reproducible. This is a **substitute source for F_MISS, not a substitute
rule**. It is consequential: the representative differs from alphabetical-only in **8 of 15**
groups, so the documented fallback would have changed which genotype survives in over half of
them. Both choices are recorded per group in `declonalization_audit.tsv`
(`is_representative` vs `rep_alphabetical`) so the difference stays measurable. Flagged for the
reviewer to accept or send back.

**3. Genericity guard on the adapter.** The first draft of `prepare_inputs.R` hardcoded the
reference/annotation/genotype paths; the repo's genericity hook flagged them. They were lifted
into `config/malay_cohort.yaml` rather than tagged as intentional hardcodes — the config is the
cohort declaration, and the lifted module takes every path as a CLI flag, so no cohort fact now
lives in code.

**Nothing skipped or blocked.** No installs attempted; `optparse` is absent but the module's
scripts hand-roll their arg parsing, as phase 0 verified.

## Proposed next (phase 2) and alternatives

Run detection for all 3 pairs × both arms under `absolute` @ 5e-4 — `Mf__Mn`, `Mf__Peninsular`,
`Mn__Peninsular` — writing `results/<run>/calls/{full,declonal}/*.tsv`. `Mf__Mn` full is already
on disk from the smoke test and will be reused rather than recomputed (same command, same
inputs). Arms differ only by `--exclude`.

Alternatives considered: (a) run the three detection rules now so phase 6's comparison is free —
rejected, the floor derivation in phase 3 reads the `absolute` calls and mixing rules into the
cached call set invites using the wrong one; the rule sweep gets its own pass on its own
scenarios. (b) Fold the arms into one run with a post-hoc sample filter — rejected, the cluster
consensus alleles are recomputed from whoever is in the cohort, so de-clonalization has to
happen *before* detection or the arm is not really de-clonalized.

## Carried forward

- The 501-vs-558 cluster-basis caveat (reviewer finding 2) is unaddressed by design and lands in
  phase 5, where the Indo 558-basis 64% ceiling and the 501 ∩ 558 overlap get reported alongside
  our coordinate reproduction.
- Peninsular is now n=28 in the de-clonalized arm, which sharpens the pre-registered risk that an
  FDR-safe global floor set by Mf leaves that cluster empty.
