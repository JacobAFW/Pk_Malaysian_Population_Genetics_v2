# MISSION — apply the agnostic introgression module to the Large-Malay cohort

## 1. Goal
Test ZOOMAL-Flow's reworked introgression module on this cohort and produce a reviewable
method result: cross-cluster introgression calls for all cluster pairs under a **floor derived
on THIS cohort**, run on both the full and de-clonalized sample sets, benchmarked **by genomic
coordinate** against the poster's 217 windows, with the detection rules compared. Deliver
Nature-Genetics-level figures and a notes doc that lets Jacob start drafting. The module lives
read-only in the Indo repo (§3); lift a copy into this project — never edit it in place.

## 2. Done when
- `results/<run>/` holds: per-pair window + sample tables (Mf/Mn/Peninsular = 3 pairs); the
  derived support-floor table (permutation + explicit FDR, this cohort — NOT Indo's N=32); a
  **coordinate-based** benchmark table vs the poster's 217 windows (CHROM + window-start, not
  WINDOW id) giving a % reproduced, for full and de-clonalized arms; a rule comparison
  (`absolute` / `relative` / `distance`) on the benchmark's measures.
- Each headline number is traceable to a named file; the clonal vs de-clonalized arms are both
  reported (do not present only one).
- Figures (vector PDF + PNG preview) in a new versioned folder; a `notes.md` summarising method,
  floor, clonality effect, benchmark, and rule choice with reasoning.

## 3. Context to read first
- `AUTONOMY_PROTOCOL.md`; this project's `MEMORY.md`; `provenance/gadi_audit.md`,
  `results/stage_b/19_diagnosis_report.md`, `docs/introgression_rewrite_brief.md` (our prior
  introgression run + the ggplot2 bug).
- Module (READ-ONLY, copy in before changing): `…/Indonesia/Pop-gen_pipeline/agnostic/` —
  `docs/introgression_analysis_spec.md` (read the **Method status** block: our "19% overlap"
  was a WINDOW-id artifact; by coordinate the density method reproduces 138/217 = 64%),
  `config/cohort.example.yaml` (introgression block), `scripts/R/introgression_*.R`, `Makefile`
  (`benchmark-malay`), `workflow/rules/05_introgression.smk`.
- Our inputs: genotype `data/processed/introgression/hmmIBD.tsv`; labels
  `data/processed/cluster_maf/cluster_labels.tsv`; GFF `data/gadi/gff/strain_A1_H.1.Icor.gff3`;
  reference `data/gadi/fasta/strain_A1_H.1.Icor.fasta`; mask `data/gadi/regions_to_mask.list`;
  clonal sets `data/raw/hmmIBD/IBD_Clonal_transmission_miss5.tsv`, `Mn_clonal_transmission.tsv`.

## 4. Suggested phases (refine in PLAN.md)
1. Copy the module into this project; assemble/verify the 5 inputs; write a Malay cohort config
   (3 clusters, `pairs:"all"`, default `absolute` @5e-4). Flag which inputs are held vs missing,
   and whether you run the lifted scripts directly or via a local Snakemake target.
2. Derive the support floor on this cohort by permutation with an explicit FDR target; record the
   floor, FDR, and its per-cluster consequence (esp. small-n Peninsular).
3. Cross-cluster detection, all 3 pairs, on BOTH full and de-clonalized sets — **headline pair
   Mf↔Mn** (the poster's intermixing finding). Write tables.
4. Coordinate-based benchmark vs the poster's 217 windows (both arms) + rule comparison
   (`absolute`/`relative`/`distance`) on Jaccard-under-drop and size-dependence.
5. Figures + `notes.md`.

## 5. Out of scope / NEEDS-JACOB
- **Focal-subgroup test = N/A (resolved from the poster).** The poster ran NO subgroup-within-
  cluster enrichment; its introgression comparison is cross-cluster pairwise — headline Mf↔Mn
  intermixing — plus an environmental association (separate run). Do NOT run the focal-enrichment
  test and do NOT invent a Malay subgroup ("Aceh" is Indo-specific). Cross-cluster only.
- NEEDS-JACOB: any floor where the FDR argument and a headline are incompatible — report as a
  finding, do not tune around it (as the spec warns for Aceh).
- Out of scope: the introgression×environment association re-run (ALC vector / tree cover /
  elevation) and its PNTD-method replication — that's a separate downstream run.
- Do not touch: the Indo `agnostic/` folder (read-only), `data/raw/` + `data/gadi/` (inputs),
  existing `figures/`, `main`.

## 6. Budgets
Phase limit 8. Max 3 threads per tool. Environment: **vvg-box** (`source vvg-box-pixi/bin/activate`)
— its R already has data.table/dplyr/sp/tibble/tidyr/tidyverse; a missing package is a blocked
action, not an install. Split the permutation floor and per-pair detection so no step exceeds ~2 h.

## 7. Output standard
Figures at Nature-Genetics level: vector PDF **plus** PNG preview, written to a NEW versioned
folder (`figures/<run>/`, `results/<run>/`); originals untouched. Palette per `scripts/_setup.R`
(viridis clusters, inferno ordinal).
