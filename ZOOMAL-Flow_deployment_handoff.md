# Handoff — run ZOOMAL-Flow on the full Large-Malay cohort

For the Large-Malay Cowork. Turn this into a Claude Code prompt that runs the
finished agnostic pop-gen pipeline (ZOOMAL-Flow) end to end on this project's
cohort. You know this dataset; this doc carries the pipeline's input contract,
which you don't. Ask back on anything underspecified rather than guessing.

## What ZOOMAL-Flow is, and its status

A data-agnostic *P. knowlesi* population-genetics pipeline (Snakemake + pixi),
refactored from Jacob's frozen V1 Indonesia pipeline into a standalone public
repo. **Complete and shipped** — all seven stages (QC, MOI/Fws, structure/
ADMIXTURE, IBD, introgression, selection, report), plus first-class clonality
handling, a runnable public example, and per-step docs.

- Repo: `https://github.com/JacobAFW/zoomal-flow` (public, MIT).
- Config-driven: you point it at a cohort with a `config.yaml` + a tidy
  `samples.tsv`. **No code edits.** "Agnostic" means no code edits, not zero
  setup — it still needs a VCF, a reference, a mask, and a per-cohort config.

**This run is an application / generalization test.** There is no ground truth
for the full Malay cohort, so it is NOT a correctness validation — correctness
is already covered by the pipeline's synthetic injected-introgression tests and
by the Malay introgression benchmark (the poster's 217 windows). The point here
is twofold: prove the finished pipeline runs clean on a full real external
cohort, and produce a genuine pop-gen result for THIS project.

## The input contract (what CC must build)

Read these two files in the repo first — they ARE the spec, with every knob
commented:

- `config/cohort.example.yaml` — the annotated template (the Indo cohort in the
  agnostic schema).
- `examples/config.yaml` — a working real config for a Malaysian ENA cohort;
  the closest analogue to what you're writing.

The pipeline needs exactly four inputs, and one tidy metadata table:

1. **A VCF** (`cohort.vcf`). Biallelic WGS SNPs. See "Which VCF" below.
2. **A reference FASTA** (`reference.fasta`) with a `.fai` — contigs are derived
   from the `.fai` minus `reference.exclude_contigs`, never hand-listed.
3. **A mask file** (`reference.mask_regions`) — you have
   `data/gadi/regions_to_mask.list`.
4. **One tidy `samples.tsv`** (`metadata.table`): one row per VCF sample, with a
   `roles` column-map naming which column is `sample_id` / `group` / `geography`
   / `country` / `host` / `date` / `case_control`. Any role left `null` just
   disables the analyses that need it (graceful degradation, logged — not a
   failure). Optional per-sample lat/long under `metadata.gis` triggers the map.

## The one true blocker to clear first — the reference

`CONTEXT-inputs.md` lists the reference as **TBD**. The pipeline cannot start
without the FASTA the Malay VCF was called against, plus its `.fai`. Before
anything else: confirm which Pk assembly `Consensus_SNPs_subset.vcf.gz` was
called against (check the VCF header contigs), locate that FASTA locally (if it
is the A1-H.1 Icor assembly, Jacob has it in the Indo project), and set
`reference.fasta` to it. Set `exclude_contigs` to match how THAT assembly names
its organelles (the Indo Icor FASTA uses `["MIT", "API"]`; the PlasmoDB build
names them by accession — check the actual `.fai`).

## The input-prep principle — DO NOT feed the drop-lists to the pipeline

ZOOMAL-Flow does only generic algorithmic operations. All cohort-specific
wrangling is done ONCE, upfront, into the clean inputs — never as pipeline
config. `data/gadi/` has several drop-lists; they split into two kinds:

- **The pipeline regenerates these itself — do NOT pre-apply them:**
  - `high_fws_exclude.txt` → Stage 2 (MOI/Fws) recomputes the high-MOI
    exclusion. Let it. Then compare its `exclude_high_moi.txt` to this file as a
    free cross-check.
  - `hmmibd_drop.txt` → Stage 4 + the clonality stage detect clonal groups and
    build the `unique` arm. Do NOT drop clonal samples upfront — the whole
    both-arms (`full` vs `unique`) comparison is the feature. Cross-check the
    detected clonal set against this file.
- **Cohort knowledge the pipeline can't derive — bake into the clean VCF /
  samples.tsv upfront:**
  - `lab_isolates.txt`, `dups_from_naming_error.txt`, `exclude_samples.txt` —
    lab strains and naming-error duplicates. Remove these samples from the VCF
    (or document them as a droplist applied in input-prep, as the Indo cohort
    does in `data/INPUT_PREP.md`) before the pipeline runs.
  - `mn_cluster_samples.txt`, the cluster CSVs — this is cluster membership;
    fold it into `samples.tsv` as the `group` column if you use reference
    labels (see below).

Controls and lane-suffix duplicates ARE handled generically in-pipeline via
`controls.exclude_patterns` and `structure.duplicate_id_pattern` — set those to
match this cohort rather than pre-filtering.

## Building samples.tsv from the six-file metadata

The Malay metadata is spread across `full_metadata.tsv`,
`Pk_clusters_metadata.csv`, `Pk_clusters_peninsular_metadata.csv`,
`full_dataset_011223.csv`, and `PK_Sabah_Sample_naming_indexes.xlsx`. The
agnostic pipeline deliberately does NOT do six-file assembly (that lived in V1
only). Collapse them upfront into ONE tidy `samples.tsv`, one row per VCF
sample, with clear columns you can map to roles. `docs/LESSONS-metadata-
wrangling.md` and `scripts/10b_attach_metadata.R` in this project are your
starting points. Verify every VCF sample ID has a metadata row before running.

## Which VCF, and the cluster-labelling choice — two calls to make

- **Which VCF:** prefer `Consensus_SNPs_subset.vcf.gz` (pre-MOI) over
  `..._no_MOI.vcf.gz`, so the pipeline does its own QC/MOI/clonality and the run
  is a true end-to-end pipeline pass (and its Fws/clonal calls become the
  cross-check above). Only use the `_no_MOI` file if you specifically want to
  match the prior analysis's post-exclusion set.
- **Cluster labelling** (`structure.cluster_labelling`): you have reference
  labels (Mn / Mf / Peninsular). Two honest options — (a) `reference`: put those
  labels in the `group` column and the pipeline uses them; fast, matches prior
  work. (b) `numbered` / `auto`: leave `group: null` and let ADMIXTURE derive
  the clusters, then cross-check against the known labels — this is the stronger
  agnosticism demonstration (pre-labelling begs the question, which is why the
  public example uses `numbered`). Recommend running derived-and-cross-checked;
  Jacob's call.

## Known chokepoints to watch (all documented in-repo)

- **ADMIXTURE determinism** is already handled (`admixture_seed` +
  `admixture_threads: 1`, per-K working dirs). Leave those as shipped.
- **hmmIBD long-path landmine:** hmmIBD has a fixed-size output-path buffer and
  SIGTRAPs past ~56 chars. Run from a SHORT working path and let the pipeline
  hand it short relative prefixes. A deep clone path is the likeliest trip.
- **Small-cohort PCA / LD-prune floor (<50 samples):** PLINK2 refuses
  `--indep-pairwise` and PCA below 50. The full Malay cohort clears this, but a
  small per-cluster subset might not — the pipeline warns and degrades rather
  than crashing (`min_samples_for_ld_prune`).
- **Both clonality arms are always produced** — expect
  `outputs/<stage>/declonalization_comparison.tsv` per relevant stage; read them.

## Guardrails (carry these into the CC prompt verbatim)

- **This is a RUN, not a repo change.** Clone `zoomal-flow` into a working area
  inside THIS project folder and run it there. No commits or pushes to
  `zoomal-flow` are needed. Do not modify the pipeline code; if something needs
  a code change, stop and report it.
- **Scratch is session-local, inside the working folder — never `/tmp`, never a
  shared path, never `mktemp`.** (A prior CC session `rm -rf`'d `/tmp/hb` and
  destroyed unrelated analysis outputs. Don't.)
- **Never delete, move, or overwrite anything you did not create this session.**
  Print any deletion target first; deletions need Jacob's explicit chat confirm.
- **Do not touch the frozen V1** Indo `Pop-gen_pipeline` or its Edy share bundle.

## What CC should report back

Step-by-step pass/fail with the exact first point anything trips; sample and SNP
counts surviving each stage; the ADMIXTURE best-K and cluster labels (+ the
cross-check against Mn/Mf/Peninsular); the clonal groups found and how the
`full` vs `unique` arms differ (+ cross-check vs `hmmibd_drop.txt`); the
introgression result; whether the report renders; and every config knob it had
to set for this cohort (so the config is reproducible). Flag anything ambiguous
rather than guessing.
