# Lessons: matching sample names across sequencing data and metadata

From the *P. knowlesi* Large Malay pop-gen project. Scope is deliberately narrow: **getting a sample in the genotype file to line up with the same sample in a metadata table.** Everything here is a failure mode we actually hit, with the concrete shape that caused it.

---

## Before writing a single join: enumerate the namespaces

- **Assume the genotype file contains more than one ID scheme.** Our `Pk.fam` holds at least three, in one column:
  - `<STUDYCODE>_DK<library>-<lane>_<flowcell>_L#` — local samples with a sequencing suffix
  - `PK_SB_DNA_<n>_DK<...>` — an entirely different prefix scheme (82 samples)
  - `ERR#######` / `SRR######` — bare public accessions (90 samples)

  One regex will not cover all three. Split the `.fam` by cohort, resolve each cohort its own way, `rbind` back — and assert the row count survives.
- **Do the same on the metadata side.** Ours came from four files keyed on `studycode`, `Sample`, `sampleid`+`group`+`subjectid`, and an ENA accession respectively. Write down the key column of each file *before* designing the join.

## Sequencing IDs are biological IDs plus machine decoration

- The seq ID is usually `<biological ID><separator><library/lane/flowcell junk>`. Ours needed `str_remove(V1, "_DK.*")`. **Strip the decoration; never try to match it.**
- Define that strip **once**, in one shared function, and use it in every script. We had the same regex chain retyped in five places; the one that forgot a step silently dropped rows rather than erroring.

## Separator, case and padding drift

- **Separators differ:** clinical metadata `studycode` is `Kud1-23` (hyphen); the fam-derived ID is `KUD1_23` (underscore). Needed `str_replace("_", "-")`.
- **Case differs:** metadata is mixed-case, fam is upper. Needed `str_to_upper()` on the metadata side.
- **PLINK rewrites your IDs.** `--vcf-idspace-to _` converts spaces in VCF sample names to underscores, and `--double-id` sets FID = IID. So a metadata field containing a space arrives in the `.fam` as an underscore. Know what your loader did to the name before you try to match it.
- **Leading zeros vanish.** The naming index stores `subjectid` as 1, 2 or 3 characters; it had to be padded to 3 before pasting to `group`. **Read every ID column as character, never as numeric** — one `read_csv` type guess turns `007` into `7` and the join dies silently.

## Some cohorts need a bridge table, not a direct join

- The `PK_SB_DNA` samples could not reach the clinical metadata at all in one hop. The path was: fam ID → strip suffix → join naming-index on `sampleid` → construct `group + zero-padded subjectid` → join clinical metadata on that. **Two hops, via a file that exists only to translate IDs.**
- Expect at least one such bridge file per project, and expect it to be an `.xlsx` someone maintained by hand.

## Public accessions: know *which* accession your data is keyed on

- ENA gives every sample several IDs: `ERR` = **run**, `ERS` = **sample**, `ERX` = experiment, `ERP`/`PRJEB` = study. Your genotype file is keyed on whichever one the FASTQs were named after.
- Ours is keyed on **run**: of the 28 accessions in the peninsular metadata, the `ERS` column matched **0** samples in the `.fam`; the `ERR` column matched **24**. Same 28 samples, same file, different column — 0% vs 86%.
- **Check the accession type explicitly, and keep a run↔sample mapping** if the metadata is published against one and the genotypes named after the other.

## Column *names* can be misaligned with column *values* — validate by shape

This is the one that bit us hardest, and it is invisible unless you look:

- `Pk_clusters_peninsular_metadata.csv` has a header row of **9** fields but data rows of **8**, because `Date collected` was written unquoted with a comma in it (`...,Date,collected,Parasitemia,...`). From that column on, **every header name is shifted one place off its values.**
- Result: the column *named* `ENA_accession_no_ES` actually contains the `ERR` (run) values, and the column *named* `ENA_accession_no_ER` is empty.
- Downstream consequence: `16v2_regen_labels.R` keys on `kp$ENA_accession_no_ER` → **zero** peninsular rows joined → the ADMIXTURE Q-column→cluster anchor had no Peninsular evidence and fell back to assign-by-elimination, while 24 genuinely matchable Peninsular IDs sat in `cleaned.fam` unused. (`10b_attach_metadata.R` selects the same column *positionally* and happens to get the right values — the two scripts disagree by accident.)
- **Lessons:**
  - Assert `length(header fields) == length(data fields)` on read. A one-field mismatch is a silent, project-wide off-by-one.
  - **Validate ID columns by shape, not by name:** does this column actually match `^ERR[0-9]+$`? If a column called `ENA_accession_no_ER` is empty and one called `Read_depth` is full of `ERS#######`, the header is lying to you.
  - Never mix positional (`select(2, 8)`) and by-name (`rename(Sample = ENA_...)`) access to the same file in the same project.

## Dirty edges: whitespace and suffixes in the key

- **Leading/trailing whitespace in headers and values.** `Pk_clusters_metadata.csv` headers are ` Code`, ` area`, ` Group` — leading spaces from Excel. We need `strip.white = TRUE` on read and `trimws()` on the key column. Add `\r` stripping for anything that has been through Windows.
- **Label columns carry inconsistent suffixes.** In the same `Group` column: `Mf-Pk`, `Mn-Pk`, and `Peninsular` — two of three values suffixed. `sub("-Pk$", "", ...)` handles it, but only because someone noticed. Tabulate the distinct values of every key and label column before using it.

## One sample ≠ one metadata row

- The clinical table is per-episode, so one `Sample` mapped to two rows. `left_join`ing it into the `.fam` inflated 755 → 756 rows and PLINK died with a duplicate-ID error.
- **Rule: `distinct(key, .keep_all = TRUE)` the metadata side before any join into a `.fam`, and assert `nrow(after) == nrow(before)` immediately after.** A `.fam` join is genotype-file surgery — row count and order are load-bearing for every downstream tool, and corruption surfaces far from the cause.

## Quantify every match — a failed join is not "missing samples"

- Print `matched / unmatched` **both directions** for every ID join, and write the unmatched IDs to a file. The match rate tells you which kind of bug you have:
  - **~0%** → wrong key type entirely (ERS vs ERR, run vs sample)
  - **~50–90%** → normalisation problem (case, separator, padding, suffix) affecting one cohort
  - **~99%** → genuinely missing metadata
- We only found the ENA bug because a diagnostic tally printed `V3: 0 known argmaxes` — a join that returns nothing looks identical to a cluster with no known members unless you count the input rows too.

## The fix that would have prevented most of this

Build **one canonical sample key table** early, and make every script join through it rather than re-deriving IDs:

```
canonical_id  fam_id  studycode  ena_run  ena_sample  cohort  source_file
```

Generate it once with all the normalisation logic in one place, check it in (or into `provenance/`), and assert on it: every `.fam` ID resolves to exactly one canonical ID, and every canonical ID appears at most once. Every later "why is n different?" question is then answered by one file.

---

### Short version

1. List every ID namespace in the genotype file *and* in each metadata file before writing a join — expect ≥3 in the `.fam` alone.
2. Seq ID = biological ID + machine suffix. Strip the suffix in one shared function.
3. Normalise case, separators and zero-padding explicitly; read all ID columns as character.
4. Know whether your data is keyed on ENA **run** or **sample** accessions — the wrong one matches 0%.
5. Validate ID columns by regex shape, not by header name; assert header field count == data field count.
6. `trimws()` keys; tabulate distinct values of key and label columns before use.
7. `distinct()` the metadata side and assert row counts before/after any join into a `.fam`.
8. Report matched/unmatched counts both ways, every time. 0% ≠ 90% ≠ 99% — each is a different bug.
9. Build one canonical sample key table and route every join through it.
