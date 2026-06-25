#!/usr/bin/env bash
# install_pk_extras.sh — layer the Pk-pop-gen-specific R packages on top of the
# vvg-box base env (see env/ENV-PLAN.md). Mirrors the Indo installer's idiom:
#   conda-available  -> pixi add (workspace)
#   Bioconductor     -> BiocManager
#   GitHub           -> remotes::install_github
#   CRAN-only        -> install.packages
#
# Prereq: Layer 0 done — vvg-box-pixi/ exists in this project root (run the Indo
# envs/install.sh first). Run from the project root. Idempotent.
#
# NOTE: this is a Mac/Claude-Code action. The Cowork Linux sandbox cannot solve
# the osx-arm64 env.

set -eo pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
cd "$PROJECT_ROOT"

if [[ ! -d "$PROJECT_ROOT/vvg-box-pixi" ]]; then
  echo "ERROR: vvg-box-pixi/ not found. Run the Layer 0 base installer first (env/ENV-PLAN.md)." >&2
  exit 1
fi

# Activate the vvg-box pixi workspace so pixi/Rscript/etc. resolve. The
# activator may dereference a symlinked vvg-box-pixi (Layer 0 is shared with
# the Indo Pop-gen install — see MEMORY.md).
ACTIVATE="$PROJECT_ROOT/vvg-box-pixi/bin/activate"
if [[ ! -f "$ACTIVATE" ]]; then
  echo "ERROR: $ACTIVATE missing. Layer 0 install incomplete." >&2
  exit 1
fi
# Tolerant sourcing — bashrc.d/96-history can be noisy on unset history files.
set +e
# shellcheck disable=SC1090
VVG_SILENT=1 source "$ACTIVATE"
set -e
command -v pixi >/dev/null || { echo "ERROR: pixi not on PATH after activation." >&2; exit 1; }

echo "==> Layer 1: pixi add Pk conda-forge R packages"
pixi workspace channel add conda-forge 2>/dev/null || true
PK_CONDA_PKGS=(
  "r-rcpp" "r-sp" "r-terra" "r-spdep" "r-geosphere"
  "r-rworldmap" "r-here" "r-car" "r-mumin"
)
echo "    pixi add ${PK_CONDA_PKGS[*]}"
pixi add "${PK_CONDA_PKGS[@]}"

echo "==> Bioconductor: LEA"
pixi run -- Rscript -e '
  if (!requireNamespace("LEA", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("LEA", ask = FALSE, update = FALSE)
  } else message("LEA already installed — skipping")
'

echo "==> GitHub: tess3r (bcm-uga/TESS3_encho_sen)"
pixi run -- Rscript -e '
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")
  if (!requireNamespace("tess3r", quietly = TRUE)) {
    remotes::install_github("bcm-uga/TESS3_encho_sen", upgrade = "never")
  } else message("tess3r already installed — skipping")
'

echo "==> CRAN: malariaAtlas, eSDM"
pixi run -- Rscript -e '
  want <- c("malariaAtlas","eSDM")
  to_install <- want[!sapply(want, requireNamespace, quietly = TRUE)]
  if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
  # eSDM fallback if archived on CRAN:
  if (!requireNamespace("eSDM", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", repos = "https://cloud.r-project.org")
    remotes::install_github("smwoodman/eSDM", upgrade = "never")
  }
'

echo "==> Layer 1 verify"
pixi run -- Rscript -e '
  pk <- c("LEA","tess3r","malariaAtlas","eSDM","terra","spdep","MuMIn","car",
          "sp","geosphere","rworldmap","here","Rcpp")
  ok <- sapply(pk, requireNamespace, quietly = TRUE)
  print(data.frame(package = pk, installed = ok))
  if (any(!ok)) { message("MISSING: ", paste(pk[!ok], collapse=", ")); quit(status=1) }
  message("All Pk extras present.")
'
echo "==> Done. Next: FEEMS env (env/ENV-PLAN.md, Layer 2)."
