#!/usr/bin/env bash
# Deterministic clean-repo scanner (vendored from JacobAFW/bio-skills scripts/scan.sh).
#
# Catches the hard-edged, pattern-detectable subset — data files, known secret
# tokens, sample sheets — that a hook can enforce automatically. The subtle
# PII/sensitivity judgement stays in the commit-guard skill; this is the part a
# machine can be sure about.
#
#   scan.sh staged    scan files staged for commit   (default)
#   scan.sh tracked   scan all tracked files          (used before a push)
#
# Exit 0 = clean. Exit 1 = findings (printed to stderr). Exit 2 = bad usage.
set -uo pipefail

mode="${1:-staged}"
case "$mode" in
  staged)  files=$(git diff --cached --name-only --diff-filter=ACM 2>/dev/null) ;;
  tracked) files=$(git ls-files 2>/dev/null) ;;
  *) echo "usage: scan.sh [staged|tracked]" >&2; exit 2 ;;
esac

findings=()

# --- filename rules ---------------------------------------------------------
data_re='\.(csv|tsv|xlsx|parquet|feather|rds|RData|rda|fastq|fq|bam|sam|cram|vcf|bcf|bed|bim|fam|ped|gff|gtf|fasta|fa|fna|h5|hdf5|npy|npz)(\.gz)?$'
secretfile_re='(^|/)(\.env($|\.)|.*\.pem$|.*\.key$|id_rsa|.*\.p12$|credentials$|\.netrc$|\.pgpass$|service-account.*\.json$)'
meta_re='(sample_?sheet|manifest|barcodes)'

while IFS= read -r f; do
  [ -z "$f" ] && continue
  printf '%s\n' "$f" | grep -Eiq "$data_re"       && findings+=("[A data]        $f")
  printf '%s\n' "$f" | grep -Eiq "$secretfile_re" && findings+=("[B secret-file] $f")
  printf '%s\n' "$f" | grep -Eiq "$meta_re"       && findings+=("[E metadata]    $f")
done <<< "$files"

# --- high-precision secret tokens in file contents --------------------------
# Only patterns with near-zero false-positive rate go here; generic keyword
# matching is left to gitleaks (below) and to the commit-guard skill.
content_re='AKIA[0-9A-Z]{16}|ghp_[0-9A-Za-z]{36}|github_pat_[0-9A-Za-z_]{20,}|-----BEGIN [A-Z ]*PRIVATE KEY-----|xox[baprs]-[0-9A-Za-z-]{10,}'
while IFS= read -r f; do
  [ -z "$f" ] && continue
  [ -f "$f" ] || continue
  line=$(grep -EnI "$content_re" "$f" 2>/dev/null | head -1)
  [ -n "$line" ] && findings+=("[B secret-content] $f:${line%%:*}")
done <<< "$files"

# --- wrap gitleaks if installed (do not reinvent secret scanning) -----------
if command -v gitleaks >/dev/null 2>&1; then
  if [ "$mode" = staged ]; then
    gitleaks protect --staged --no-banner >/dev/null 2>&1 || findings+=("[B gitleaks] staged changes flagged by gitleaks")
  else
    gitleaks detect --no-banner >/dev/null 2>&1 || findings+=("[B gitleaks] history flagged by gitleaks")
  fi
fi

if [ "${#findings[@]}" -gt 0 ]; then
  echo "clean-repo scan: ${#findings[@]} issue(s) in the $mode set" >&2
  printf '  %s\n' "${findings[@]}" >&2
  echo "Fix before continuing: gitignore + 'git rm --cached' for data/metadata; rotate-then-scrub for secrets." >&2
  exit 1
fi
exit 0
