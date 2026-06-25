# Large Malay *P. knowlesi* population genetics

Downstream population-genetics analysis code for *Plasmodium knowlesi* whole-genome
sequencing from Malaysia (Sabah + Peninsular).

## What this repo is

Analysis scripts that take a called VCF through population structure, relatedness,
selection, introgression, and migration-surface analysis for *P. knowlesi*. The
upstream calling pipeline runs on HPC (NCI Gadi) and is out of scope here; this repo
is the local, downstream analysis only. Stack: R / R Markdown + PBS job scripts,
reproduced with a pixi environment; FEEMS (Python) for the migration surface.

## What's included — and what's deliberately not

This repository contains **code only**. By design it does **not** include:

- raw or processed **data** (the VCFs, reference genome, GFF/BED, genotype tables)
- **sample sheets, manifests, or metadata** that link samples to individuals
  (e.g. `PK_Sabah_Sample_naming_indexes.xlsx`, cluster metadata, GPS coordinates)
- any **identifying or sensitive** information (sample-collection coordinates;
  the unpublished conference poster)
- credentials, tokens, or environment files; vendored binaries and the ~6 GB
  pixi/conda environment trees

Data lives on the institutional HPC (NCI Gadi) and is not in version control. The
scripts expect inputs under `data/raw/` at the paths described below.

> No secrets or credentials were found in the included files. One redaction is
> still required before publishing — see **Before publishing**.

## Reproducing the analysis

1. **Environment.** Core tools (bcftools, samtools, plink/plink2, R 4.5, quarto)
   from `env/pixi.toml.indo-base` + `env/pixi.lock.indo-base` via `pixi install`.
   Then R extras via `env/install_pk_extras.sh`. FEEMS is a separate micromamba
   env; apply `env/feems_numpy2_patch.sh` (NumPy-2 fix). FEEMS requires
   **segregating sites only** — drop monomorphic sites before `fit()`.
2. **Inputs.** Place the consensus SNP VCF (bgzipped + indexed), the
   *P. knowlesi* A1-H.1 reference FASTA (+ `.fai`/`.dict`) and GFF, and sample
   metadata (with coordinates for FEEMS) under `data/raw/`.
3. **Run.** Work through `scripts/01_pop_gen.Rmd` → `02_spatial_introgression.Rmd`
   → `03_eems_prep.Rmd`; HPC stages (hmmIBD, rehh, ADMIXTURE, TESS3, ordination,
   introgression) are in `scripts/gadi/`.

## Structure

```
scripts/        # analysis code: top-level Rmd workflow + scripts/gadi/ (R + PBS)
env/            # pixi manifest + lockfile, R-extras and FEEMS patch scripts
data/           # (empty) expected input location — not tracked
results/  figures/  outputs/   # (empty) output targets — not tracked
```

## Before publishing

- **Redact the HPC username `jw1542`** (`/home/588/jw1542/...`) in
  `scripts/gadi/ibd_and_selection.Rmd`, `scripts/gadi/Analyses.Rmd`, and
  `scripts/gadi/selection_cluster_spec.pbs` — category-D identifying info.
- **Generalise absolute paths** (`/g/data/pq84/...`, `/Users/.../Pk_Pipeline/...`)
  in `scripts/gadi/` to relative `data/raw/...` paths so the repo runs from a
  fresh clone.

## License

<add a license before publishing — e.g. MIT for code>
